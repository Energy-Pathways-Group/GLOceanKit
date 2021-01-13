classdef Boussinesq3DConstantStratification < handle
    %3D Boussinesq model with constant stratification solved in wave-vortex
    %space
    
    properties
        x, y, z
        k, l, j
        Nx, Ny, Nz
        Lx, Ly, Lz
        f0, N0
        
        dctScratch, dstScratch;
        F,G
        ApU, ApV, ApN
        AmU, AmV, AmN
        A0U, A0V, A0N
        
        UAp, UAm, UA0
        VAp, VAm, VA0
        WAp, WAm
        NAp, NAm, NA0
    end
    
    properties (Dependent)
        % These convert the coefficients to their depth integrated energies
        Apm_HKE_factor
        Apm_VKE_factor
        Apm_PE_factor
        Apm_TE_factor

        A0_HKE_factor
        A0_PE_factor
        A0_TE_factor
    end
    
    properties (Constant)
        g = 9.81;
    end
    
    methods
        function self = Boussinesq3DConstantStratification(dims, n, latitude, N0, rho0)
            % rho0 is optional.
            if length(dims) ~=3 || length(n) ~= 3
                error('The dims and n variables must be of length 3. You need to specify x,y,z');
            end
            
            if mod(log2(n(3)),1) == 0
                error('You are implicitly asking for periodic boundary conditions in the vertical. This is not supported.');
            elseif mod(log2(n(3)-1),1) == 0
                self.Nz = n(3);
            else
                error('The vertical dimension must have 2^n or (2^n)+1 points. This is an artificial restriction.');
            end
                        
            self.Lx = dims(1);
            self.Ly = dims(2);
            self.Lz = dims(3);
            
            self.Nx = n(1);
            self.Ny = n(2);
            self.Nz = n(3);

            dx = self.Lx/self.Nx;
            dy = self.Ly/self.Ny;
            dz = self.Lz/(self.Nz-1);
            
            self.x = dx*(0:self.Nx-1)'; % periodic basis
            self.y = dy*(0:self.Ny-1)'; % periodic basis
            self.z = dz*(0:(self.Nz-1))' - self.Lz; % Cosine basis for DCT-I and DST-I
                        
            dk = 1/self.Lx;          % fourier frequency
            self.k = 2*pi*([0:ceil(self.Nx/2)-1 -floor(self.Nx/2):-1]*dk)';
            dl = 1/self.Ly;          % fourier frequency
            self.l = 2*pi*([0:ceil(self.Ny/2)-1 -floor(self.Ny/2):-1]*dl)';
            self.j = (0:(self.Nz-2))';
            
            self.N0 = N0;
            self.f0 = 2 * 7.2921E-5 * sin( latitude*pi/180 );
            
            self.BuildTransformationMatrices();
            
            % Preallocate this array for a faster dct
            self.dctScratch = zeros(self.Nx,self.Ny,2*(self.Nz-1));
            self.dstScratch = complex(zeros(self.Nx,self.Ny,2*(self.Nz-1)));
        end
        
        function self = BuildTransformationMatrices(self)
            % Build wavenumbers
            [K,L,J] = ndgrid(self.k,self.l,self.j);
            alpha = atan2(L,K);
            K2 = K.*K + L.*L;
            Kh = sqrt(K2);      % Total horizontal wavenumber
            M = J*pi/self.Lz;        % Vertical wavenumber
            
            N = self.N0;
            f = self.f0;
            
            g_ = 9.81;
            h = (1/g_)*(N*N-f*f)./(M.*M+K2);
            h(:,:,1) = 1; % prevent divide by zero
            
            omega = sqrt(g_*h.*K2 + f*f);
            fOmega = f./omega;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Normalization for the vertical modes
            % This comes from equations B12 in the manuscript.
            signNorm = -2*(mod(J,2) == 1)+1;
            self.F = signNorm .* (h.*M)*sqrt(2*g_/(self.Lz*(N*N-f*f)));
            self.G = signNorm .* sqrt(2*g_/(self.Lz*(N*N-f*f)));
            self.F(:,:,1) = 1; % j=0 mode
            self.G(:,:,1) = 1; % j=0 mode
            
            MakeHermitian = @(f) InternalWaveModel.MakeHermitian(f);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Transform matrices (U,V,N) -> (Ap,Am,A0)
            % This comes from equations B13 and B14 in the manuscript
            % or equation 5.5 without the factor of h.
            self.ApU = (1/2)*(cos(alpha)+sqrt(-1)*fOmega.*sin(alpha));
            self.ApV = (1/2)*(sin(alpha)-sqrt(-1)*fOmega.*cos(alpha));
            self.ApN = -g_*Kh./(2*omega);
            
            self.AmU = (1/2)*(cos(alpha)-sqrt(-1)*fOmega.*sin(alpha));
            self.AmV = (1/2)*(sin(alpha)+sqrt(-1)*fOmega.*cos(alpha));
            self.AmN = g_*Kh./(2*omega);
            
            % There are no k^2+l^2>0, j=0 wave solutions. Only the inertial
            % solution exists at k=l=j=0.
            self.ApU(:,:,1) = 0;
            self.ApV(:,:,1) = 0;
            self.ApN(:,:,1) = 0;
            
            self.AmU(:,:,1) = 0;
            self.AmV(:,:,1) = 0;
            self.AmN(:,:,1) = 0;
            
            % Now set the inertial stuff (this is just a limit of above)
            self.ApU(1,1,:) = 1/2;
            self.ApV(1,1,:) = -sqrt(-1)/2;
            self.AmU(1,1,:) = 1/2;
            self.AmV(1,1,:) = sqrt(-1)/2;
            
            % Equation B14
            self.A0U = sqrt(-1)*h.*(fOmega./omega) .* L;
            self.A0V = -sqrt(-1)*h.*(fOmega./omega) .* K;
            self.A0N = fOmega.^2;
            
            % k > 0, l > 0, j=0; Equation B11 in the manuscript
            self.A0U(:,:,1) =  sqrt(-1)*(f/g_)*L(:,:,1)./K2(:,:,1); % Note the divide by zero at k=l=0
            self.A0V(:,:,1) = -sqrt(-1)*(f/g_)*K(:,:,1)./K2(:,:,1);
            self.A0N(:,:,1) = 0;
            
            % There are no k=l=0, j>=0 geostrophic solutions
            self.A0U(1,1,:) = 0;
            self.A0V(1,1,:) = 0;
            self.A0N(1,1,:) = 0;
            
            % Finally, we need to take care of the extra factor of 2 that
            % comes out of the discrete cosine transform
            
            % Now make the Hermitian conjugate match.
            self.ApU = MakeHermitian(self.ApU);
            self.ApV = MakeHermitian(self.ApV);
            self.ApN = MakeHermitian(self.ApN);
            
            self.AmU = MakeHermitian(self.AmU);
            self.AmV = MakeHermitian(self.AmV);
            self.AmN = MakeHermitian(self.AmN);
            
            self.A0U = MakeHermitian(self.A0U);
            self.A0V = MakeHermitian(self.A0V);
            self.A0N = MakeHermitian(self.A0N);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Transform matrices (Ap,Am,A0) -> (U,V,W,N)
            % These can be pulled from equation C4 in the manuscript
            self.UAp = (cos(alpha)-sqrt(-1)*fOmega.*sin(alpha));
            self.UAm = (cos(alpha)+sqrt(-1)*fOmega.*sin(alpha));
            self.UA0 = -sqrt(-1)*(g_/f)*L;
            
            self.VAp = (sin(alpha)+sqrt(-1)*fOmega.*cos(alpha));
            self.VAm = (sin(alpha)-sqrt(-1)*fOmega.*cos(alpha));
            self.VA0 = sqrt(-1)*(g_/f)*K;
                
            self.WAp = -sqrt(-1)*Kh.*h;
            self.WAm = -sqrt(-1)*Kh.*h;
            
            self.NAp = -Kh.*h./omega;
            self.NAm = Kh.*h./omega;
            self.NA0 = ones(size(Kh));
            
            % No buoyancy anomaly for j=0 geostrophic solutions
            self.NA0(:,:,1) = 0;
            
            % There are no k^2+l^2>0, j=0 wave solutions. 
            self.UAp(:,:,1) = 0;
            self.VAp(:,:,1) = 0;
            self.NAp(:,:,1) = 0;
            
            self.UAm(:,:,1) = 0;
            self.VAm(:,:,1) = 0;
            self.NAm(:,:,1) = 0;
            
            % Only the inertial solution exists at k=l=j=0 as a negative
            % wave.
            self.UAp(1,1,:) = 1;
            self.VAp(1,1,:) = sqrt(-1);
            self.UAm(1,1,:) = 1;
            self.VAm(1,1,:) = -sqrt(-1);
            
            % Now make the Hermitian conjugate match.
            self.UAp = MakeHermitian(self.UAp);
            self.UAm = MakeHermitian(self.UAm);
            self.UA0 = MakeHermitian(self.UA0);
            
            self.VAp = MakeHermitian(self.VAp);
            self.VAm = MakeHermitian(self.VAm);
            self.VA0 = MakeHermitian(self.VA0);
            
            self.WAp = MakeHermitian(self.WAp);
            self.WAm = MakeHermitian(self.WAm);
            
            self.NAp = MakeHermitian(self.NAp);
            self.NAm = MakeHermitian(self.NAm);
            self.NA0 = MakeHermitian(self.NA0);
        end
        
        function value = get.Apm_TE_factor(self)
            [K,L,J] = ndgrid(self.k,self.l,self.j);
            K2 = K.*K + L.*L;
            M = J*pi/self.Lz;
            h = (1/self.g)*(self.N0*self.N0-self.f0*self.f0)./(M.*M+K2);
            
            value = h; % factor of 2 larger than in the manuscript
            value(:,:,1) = self.Lz/4; % factor of 4 smaller, to account for the j=0 scaling of the DCT-I
        end
        
        function value = get.A0_HKE_factor(self)
            [K,L,J] = ndgrid(self.k,self.l,self.j);
            K2 = K.*K + L.*L;
            M = J*pi/self.Lz;
            h = (1/self.g)*(self.N0*self.N0-self.f0*self.f0)./(M.*M+K2);
            h(:,:,1) = 1; % prevent divide by zero 
            %omega = sqrt(self.g*h.*K2 + self.f0*self.f0);
            
            % This comes from equation (3.10) in the manuscript, but using
            % the relation from equation A2b
%             value = (self.g/(self.f0*self.f0)) * (omega.*omega - self.f0*self.f0) .* (self.N0*self.N0 - omega.*omega) / (2 * (self.N0*self.N0 - self.f0*self.f0) );
            value = (self.g^3/(self.f0*self.f0)) * K2.*h.*h.*M.*M / (2 * (self.N0*self.N0 - self.f0*self.f0) ); % factor of 2 larger than in the manuscript
            value(:,:,1) = (self.g^2/(self.f0*self.f0)) * K2(:,:,1) * self.Lz/8; % factor of 4 smaller, to account for the j=0 scaling of the DCT-I
        end
        function value = get.A0_PE_factor(self)
            value = self.g*self.N0*self.N0/(self.N0*self.N0-self.f0*self.f0)/2; % factor of 2 larger than in the manuscript
        end
        function value = get.A0_TE_factor(self)
            value = self.A0_HKE_factor + self.A0_PE_factor;
        end
                
        function [Ap,Am,A0] = Project(self,U,V,N)
            Ubar = self.TransformFromSpatialDomainWithF(U);
            Vbar = self.TransformFromSpatialDomainWithF(V);
            Nbar = self.TransformFromSpatialDomainWithG(N);
            
            Ap = self.ApU.*Ubar + self.ApV.*Vbar + self.ApN.*Nbar;
            Am = self.AmU.*Ubar + self.AmV.*Vbar + self.AmN.*Nbar;
            A0 = self.A0U.*Ubar + self.A0V.*Vbar + self.A0N.*Nbar;
        end
        
        function [U,V,W,N] = VelocityField(self,Ap,Am,A0)
            Ubar = self.UAp.*Ap + self.UAm.*Am + self.UA0.*A0;
            Vbar = self.VAp.*Ap + self.VAm.*Am + self.VA0.*A0;
            Wbar = self.WAp.*Ap + self.WAm.*Am;
            Nbar = self.NAp.*Ap + self.NAm.*Am + self.NA0.*A0;
            
            U = self.TransformToSpatialDomainWithF(Ubar);
            V = self.TransformToSpatialDomainWithF(Vbar);
            W = self.TransformToSpatialDomainWithG(Wbar);
            N = self.TransformToSpatialDomainWithG(Nbar);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Transformations to and from the spatial domain
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [C11,C21,C31,C12,C22,C32,C13,C23,C33] = Validate(self)
            % This is S*S^{-1} and therefore returns the values in
            % wave-vortex space. So, C11 represents Ap and should be 1s
            % where we expected Ap solutions to exist.
            C11 = self.ApU.*self.UAp + self.ApV.*self.VAp + self.ApN.*self.NAp;
            C21 = self.AmU.*self.UAp + self.AmV.*self.VAp + self.AmN.*self.NAp;
            C31 = self.A0U.*self.UAp + self.A0V.*self.VAp + self.A0N.*self.NAp;
            
            C12 = self.ApU.*self.UAm + self.ApV.*self.VAm + self.ApN.*self.NAm;
            C22 = self.AmU.*self.UAm + self.AmV.*self.VAm + self.AmN.*self.NAm;
            C32 = self.A0U.*self.UAm + self.A0V.*self.VAm + self.A0N.*self.NAm;
            
            C13 = self.ApU.*self.UA0 + self.ApV.*self.VA0 + self.ApN.*self.NA0;
            C23 = self.AmU.*self.UA0 + self.AmV.*self.VA0 + self.AmN.*self.NA0;
            C33 = self.A0U.*self.UA0 + self.A0V.*self.VA0 + self.A0N.*self.NA0;
            
            % Maybe check that max(abs(C12(:))) is very small.
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Transformations to and from the spatial domain
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        
        function u = TransformToSpatialDomainWithF(self, u_bar)
            u = self.Nx*self.Ny*ifft(ifft(u_bar.*self.F,self.Nx,1),self.Ny,2,'symmetric');
            self.dctScratch = cat(3, 0.5*u(:,:,1:self.Nz-1), zeros(self.Nx,self.Ny), 0.5*u(:,:,(self.Nz-1):-1:2));
            u = fft(self.dctScratch,2*(self.Nz-1),3);

            % should not have to call real, but for some reason, with enough
            % points, it starts generating some small imaginary component.
            u = real(u(:,:,1:self.Nz));
        end
        
        function w = TransformToSpatialDomainWithG(self, w_bar )
            % Here we use what I call the 'Fourier series' definition of the ifft, so
            % that the coefficients in frequency space have the same units in time.
            w = self.Nx*self.Ny*ifft(ifft(w_bar.* self.G,self.Nx,1),self.Ny,2,'symmetric');
            
            % Re-order to convert to an fast cosine transform
            self.dstScratch = 0.5*sqrt(-1)*cat(3, zeros(self.Nx,self.Ny), w(:,:,2:(self.Nz-1)), zeros(self.Nx,self.Ny), -w(:,:,(self.Nz-1):-1:2));
            
            w = fft( self.dstScratch,2*(self.Nz-1),3);
            % should not have to call real, but for some reason, with enough
            % points, it starts generating some small imaginary component.
            w = real(w(:,:,1:self.Nz)); 
        end
        
        function u_bar = TransformFromSpatialDomainWithF(self, u)
            % df = 1/(2*(Nz-1)*dz)
            % nyquist = (Nz-1)*df
            self.dctScratch = ifft(cat(3,u,u(:,:,(self.Nz-1):-1:2)),2*(self.Nz-1),3);
            u_bar = 2*real(self.dctScratch(:,:,1:(self.Nz-1))); % include barotropic mode, but leave off the Nyquist.
            u_bar = fft(fft(u_bar,self.Nx,1),self.Ny,2)./self.Nx/self.Ny;
            u_bar = u_bar./self.F;
        end
        
        function w_bar = TransformFromSpatialDomainWithG(self, w)
            % df = 1/(2*(Nz-1)*dz)
            % nyquist = (Nz-2)*df
            self.dstScratch = ifft(cat(3,w,-w(:,:,(self.Nz-1):-1:2)),2*(self.Nz-1),3);
            w_bar = 2*imag(self.dstScratch(:,:,1:(self.Nz-1)));
            w_bar = fft(fft(w_bar,self.Nx,1),self.Ny,2)/self.Nx/self.Ny;
            w_bar = w_bar./self.G;
        end
    end
end

