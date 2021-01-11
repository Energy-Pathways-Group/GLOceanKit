classdef Boussinesq3DConstantStratification
    %3D Boussinesq model with constant stratification solved in wave-vortex
    %space
    
    properties
        x, y, z
        Nx, Ny, Nz, nz
        dctScratch, dstScratch;
        F,G
        ApU, ApV, ApN
        AmU, AmV, AmN
        A0U, A0V, A0N
    end
    
    methods
        function self = Boussinesq3DConstantStratification(dims, n, latitude, N0, rho0)
            % rho0 is optional.
            if length(dims) ~=3 || length(n) ~= 3
                error('The dims and n variables must be of length 3. You need to specify x,y,z');
            end
            
            if mod(log2(n(3)),1) == 0
                error('You are implicitly asking for periodic boundary conditions in the vertical. This is not supported.');
            elseif mod(log2(n(3)-1),1) == 0 % user wants the surface point
                nz = n(3)-1; % internally we proceed as if there are n-1 points
            else
                error('The vertical dimension must have 2^n or (2^n)+1 points. This is an artificial restriction.');
            end
            
            % Construct the vertical dimension
            Lz = dims(3);
            Nz = n(3);
            dz = Lz/nz;
            z = dz*(0:Nz-1)' - Lz;
            
            Lx = dims(1);
            Ly = dims(2);
            
            self.Nx = n(1);
            self.Ny = n(2);
            self.Nz = n(3);
            self.nz = nz;
            
            % Preallocate this array for a faster dct
            self.dctScratch = zeros(self.Nx,self.Ny,2*nz);
            self.dstScratch = complex(zeros(self.Nx,self.Ny,2*nz));
            
            dx = Lx/self.Nx;
            dy = Ly/self.Ny;
            
            self.x = dx*(0:self.Nx-1)'; % periodic basis
            self.y = dy*(0:self.Ny-1)'; % periodic basis
            self.z = z; % cosine basis (not your usual dct basis, however)
                        
            dk = 1/Lx;          % fourier frequency
            k = 2*pi*([0:ceil(self.Nx/2)-1 -floor(self.Nx/2):-1]*dk)';
            dl = 1/Ly;          % fourier frequency
            l = 2*pi*([0:ceil(self.Ny/2)-1 -floor(self.Ny/2):-1]*dl)';
            j = (1:(nz-1))';
            
            [K,L,J] = ndgrid(k,l,j);
            alpha = atan2(L,K);
            f0 = 2 * 7.2921E-5 * sin( latitude*pi/180 );
            K2 = K.*K + L.*L;
            Kh = sqrt(K2);
            
            g = 9.81;
            M = J*pi/Lz;        % Vertical wavenumber
            h = (1/g)*(N0*N0-f0*f0)./(M.*M+K2);
            omega = sqrt(g*h.*K2 + f0*f0);
            fOmega = f0./omega;
            
            signNorm = -2*(mod(J,2) == 1)+1;
            self.F = signNorm .* (h.*M)*sqrt(2*g/(Lz*(N0*N0-f0*f0)));
            self.G = signNorm .* sqrt(2*g/(Lz*(N0*N0-f0*f0)));
            
            self.ApU = InternalWaveModel.MakeHermitian((1/2)*(cos(alpha)+sqrt(-1)*fOmega.*sin(alpha)));
            self.ApV = InternalWaveModel.MakeHermitian((1/2)*(sin(alpha)-sqrt(-1)*fOmega.*cos(alpha)));
            self.ApN = InternalWaveModel.MakeHermitian(-g*Kh./(2*omega));
            
            self.AmU = InternalWaveModel.MakeHermitian((1/2)*(cos(alpha)-sqrt(-1)*fOmega.*sin(alpha)));
            self.AmV = InternalWaveModel.MakeHermitian((1/2)*(sin(alpha)+sqrt(-1)*fOmega.*cos(alpha)));
            self.AmN = InternalWaveModel.MakeHermitian(g*Kh./(2*omega));
            
            self.A0U = InternalWaveModel.MakeHermitian(sqrt(-1)*h.*(fOmega./omega) .* L);
            self.A0V = InternalWaveModel.MakeHermitian(-sqrt(-1)*h.*(fOmega./omega) .* K);
            self.A0N = InternalWaveModel.MakeHermitian(fOmega.^2);
        end
        
        function [Ap,Am,A0] = Project(self,U,V,N)
            Ubar = self.TransformFromSpatialDomainWithF(U);
            Vbar = self.TransformFromSpatialDomainWithF(V);
            Nbar = self.TransformFromSpatialDomainWithG(N);
            
            Ap = self.ApU.*Ubar + self.ApV.*Vbar + self.ApN.*Nbar;
            Am = self.AmU.*Ubar + self.AmV.*Vbar + self.AmN.*Nbar;
            A0 = self.A0U.*Ubar + self.A0V.*Vbar + self.A0N.*Nbar;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computes the phase information given the amplitudes (internal)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function u = TransformToSpatialDomainWithBarotropicFMode(self, u_bar)
            u = self.Nx*self.Ny*ifft(ifft(u_bar,self.Nx,1),self.Ny,2,'symmetric');
            u = repmat(u,[1 1 self.Nz]);
        end
        
        function u = TransformToSpatialDomainWithF(self, u_bar)
            % Note that the extra multiplication by self.F could be
            % precomputed---but we keep it for simplicity, at the cost of
            % speed.
            
            % Here we use what I call the 'Fourier series' definition of the ifft, so
            % that the coefficients in frequency space have the same units in time.
            u = self.Nx*self.Ny*ifft(ifft(u_bar.*self.F,self.Nx,1),self.Ny,2,'symmetric');
            
            % Re-order to convert to an fast cosine transform
            % There's a lot of time spent here, and below taking the real
            % part.
            self.dctScratch = cat(3, zeros(self.Nx,self.Ny), 0.5*u(:,:,1:self.nz-1), zeros(self.Nx,self.Ny), 0.5*u(:,:,self.nz-1:-1:1));
            %             self.dctScratch(:,:,2:self.nz) =  0.5*u(:,:,1:self.nz-1); % length= nz-1
            %             self.dctScratch(:,:,self.nz+1) =  u(:,:,self.nz); % length 1
            %             self.dctScratch(:,:,(self.nz+2):(2*self.nz)) = 0.5*u(:,:,self.nz-1:-1:1); % length nz-1
            
            u = fft(self.dctScratch,2*self.nz,3);
            if self.performSanityChecks == 1
                ratio = max(max(max(abs(imag(u)))))/max(max(max(abs(real(u)))));
                if ratio > 1e-6
                    fprintf('WARNING: The inverse cosine transform reports an unreasonably large imaginary part, %.2g.\n',ratio);
                end
            end
            % should not have to call real, but for some reason, with enough
            % points, it starts generating some small imaginary component.
            u = real(u(:,:,1:self.Nz)); % Here we use Nz (not nz) because the user may want the end point.
        end
        
        function w = TransformToSpatialDomainWithG(self, w_bar )
            % Here we use what I call the 'Fourier series' definition of the ifft, so
            % that the coefficients in frequency space have the same units in time.
            w = self.Nx*self.Ny*ifft(ifft(w_bar.* self.G,self.Nx,1),self.Ny,2,'symmetric');
            
            % Re-order to convert to an fast cosine transform
            self.dstScratch = sqrt(-1)*cat(3, zeros(self.Nx,self.Ny), 0.5*w(:,:,1:self.nz-1), zeros(self.Nx,self.Ny), -0.5*w(:,:,self.nz-1:-1:1));
            
            w = fft( self.dstScratch,2*self.nz,3);
            if self.performSanityChecks == 1
                ratio = max(max(max(abs(imag(w)))))/max(max(max(abs(real(w)))));
                if ratio > 1e-6
                    fprintf('WARNING: The inverse sine transform reports an unreasonably large imaginary part, %.2g.\n',ratio);
                end
            end
            % should not have to call real, but for some reason, with enough
            % points, it starts generating some small imaginary component.
            w = real(w(:,:,1:self.Nz)); % Here we use Nz (not nz) because the user may want the end point.
        end
        
        function u_bar = TransformFromSpatialDomainWithBarotropicFMode(self, u)
            % Consistent with the DCT-I, the end points only have half the
            % width of the other points.
            u(:,:,1) = 0.5*u(:,:,1);
            u(:,:,end) = 0.5*u(:,:,end);
            u_bar = fft(fft(sum(u,3)/(self.Nz-1),self.Nx,1),self.Ny,2)/self.Nx/self.Ny;
        end
        
        function u_bar = TransformFromSpatialDomainWithF(self, u)
            % df = 1/(2*(Nz-1)*dz)
            % nyquist = (Nz-1)*df
            self.dctScratch = ifft(cat(3,u,u(:,:,self.nz:-1:2)),2*self.nz,3);
            u_bar = 2*real(self.dctScratch(:,:,2:self.nz)); % we *ignore* the barotropic mode, starting at 2, instead of 1 in the transform
            u_bar = fft(fft(u_bar,self.Nx,1),self.Ny,2)./self.F/self.Nx/self.Ny;
        end
        
        function w_bar = TransformFromSpatialDomainWithG(self, w)
            % df = 1/(2*(Nz-1)*dz)
            % nyquist = (Nz-2)*df
            self.dstScratch = ifft(cat(3,w,-w(:,:,self.nz:-1:2)),2*self.nz,3);
            w_bar = 2*imag(self.dstScratch(:,:,2:self.nz));
            w_bar = fft(fft(w_bar,self.Nx,1),self.Ny,2)./self.G/self.Nx/self.Ny;
        end
    end
end

