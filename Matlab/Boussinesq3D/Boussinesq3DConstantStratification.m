classdef Boussinesq3DConstantStratification < handle
    %3D Boussinesq model with constant stratification solved in wave-vortex
    %space
    
    properties
        x, y, z
        k, l, j
        Nx, Ny, Nz
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
                        
            Lx = dims(1);
            Ly = dims(2);
            Lz = dims(3);
            
            self.Nx = n(1);
            self.Ny = n(2);
            self.Nz = n(3);

            dx = Lx/self.Nx;
            dy = Ly/self.Ny;
            dz = Lz/(self.Nz-1);
            
            self.x = dx*(0:self.Nx-1)'; % periodic basis
            self.y = dy*(0:self.Ny-1)'; % periodic basis
            self.z = dz*(0:(self.Nz-1))' - Lz; % Cosine basis for DCT-I and DST-I
                        
            dk = 1/Lx;          % fourier frequency
            self.k = 2*pi*([0:ceil(self.Nx/2)-1 -floor(self.Nx/2):-1]*dk)';
            dl = 1/Ly;          % fourier frequency
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
                Lz = max(self.z)-min(self.z);
                M = J*pi/Lz;        % Vertical wavenumber
                
                N = self.N0;
                f = self.f0;
                
                g = 9.81;
                h = (1/g)*(N*N-f*f)./(M.*M+K2);
                h(:,:,1) = 1; % prevent divide by zero 
                
                omega = sqrt(g*h.*K2 + f*f);
                fOmega = f./omega;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Normalization for the vertical modes
                % This comes from equations B12 in the manuscript.
                signNorm = -2*(mod(J,2) == 1)+1;
                self.F = signNorm .* (h.*M)*sqrt(2*g/(Lz*(N*N-f*f)));
                self.G = signNorm .* sqrt(2*g/(Lz*(N*N-f*f)));
                self.F(:,:,1) = 1; % j=0 mode
                self.G(:,:,1) = 1; % j=0 mode
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Transform matrices (U,V,N) -> (Ap,Am,A0)
                % This comes from equations B13 and B14 in the manuscript
                self.ApU = InternalWaveModel.MakeHermitian((1/2)*(cos(alpha)+sqrt(-1)*fOmega.*sin(alpha)));
                self.ApV = InternalWaveModel.MakeHermitian((1/2)*(sin(alpha)-sqrt(-1)*fOmega.*cos(alpha)));
                self.ApN = InternalWaveModel.MakeHermitian(-g*Kh./(2*omega));
                
                self.AmU = InternalWaveModel.MakeHermitian((1/2)*(cos(alpha)-sqrt(-1)*fOmega.*sin(alpha)));
                self.AmV = InternalWaveModel.MakeHermitian((1/2)*(sin(alpha)+sqrt(-1)*fOmega.*cos(alpha)));
                self.AmN = InternalWaveModel.MakeHermitian(g*Kh./(2*omega));
                
                self.A0U = InternalWaveModel.MakeHermitian(sqrt(-1)*h.*(fOmega./omega) .* L);
                self.A0V = InternalWaveModel.MakeHermitian(-sqrt(-1)*h.*(fOmega./omega) .* K);
                self.A0N = InternalWaveModel.MakeHermitian(fOmega.^2);
                
                % k=l=0, j>=0 has only an inertial solution, which is
                % correctly created using the internal wave relations.
                % However, this means we need zero the j=0 geostrophic mode
                self.A0U(1,1,:) = 0;
                self.A0V(1,1,:) = 0;
                self.A0N(1,1,:) = 0;
                
                % k > 0, l > 0, j=0; Equation B11 in the manuscript
                self.A0U(:,:,1) =  sqrt(-1)*(f/g)*L(:,:,1)./K2(:,:,1);
                self.A0V(:,:,1) = -sqrt(-1)*(f/g)*K(:,:,1)./K2(:,:,1);
                self.A0N(:,:,1) = 0;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Transform matrices (Ap,Am,A0) -> (U,V,W,N)
                % These can be pulled from equation C4 in the manuscript
                self.UAp = InternalWaveModel.MakeHermitian((cos(alpha)-sqrt(-1)*fOmega.*sin(alpha)));
                self.UAm = InternalWaveModel.MakeHermitian((cos(alpha)+sqrt(-1)*fOmega.*sin(alpha)));
                self.UA0 = InternalWaveModel.MakeHermitian(-sqrt(-1)*(g/f)*L);
                
                self.VAp = InternalWaveModel.MakeHermitian((sin(alpha)+sqrt(-1)*fOmega.*cos(alpha)));
                self.VAm = InternalWaveModel.MakeHermitian((sin(alpha)-sqrt(-1)*fOmega.*cos(alpha)));
                self.VA0 = InternalWaveModel.MakeHermitian(sqrt(-1)*(g/f)*K);
                
                self.WAp = InternalWaveModel.MakeHermitian(-sqrt(-1)*Kh.*h);
                self.WAm = InternalWaveModel.MakeHermitian(sqrt(-1)*Kh.*h);
                
                self.NAp = InternalWaveModel.MakeHermitian(-Kh.*h./omega);
                self.NAm = InternalWaveModel.MakeHermitian(Kh.*h./omega);
                self.NA0 = ones(size(Kh));
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
            W = self.TransformToSpatialDomainWithF(Wbar);
            N = self.TransformToSpatialDomainWithF(Nbar);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computes the phase information given the amplitudes (internal)
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

