classdef WaveVortexModelConstantStratification < WaveVortexModel
    %3D Boussinesq model with constant stratification solved in wave-vortex
    %space
    
    properties
        N0   
        realScratch, complexScratch; % of size Nx x Ny x (2*Nz-1)
        F,G
        
        h
        Apm_TE_factor
        A0_HKE_factor
        A0_PE_factor
        A0_TE_factor
    end
        
    methods
        function self = WaveVortexModelConstantStratification(dims, n, latitude, N0, rho0)
            % rho0 is optional.
            if length(dims) ~=3 || length(n) ~= 3
                error('The dims and n variables must be of length 3. You need to specify x,y,z');
            end
            
            if mod(log2(n(3)),1) == 0
                error('You are implicitly asking for periodic boundary conditions in the vertical. This is not supported.');
            elseif mod(log2(n(3)-1),1) == 0
                Nz = n(3);
            else
                error('The vertical dimension must have 2^n or (2^n)+1 points. This is an artificial restriction.');
            end
                        

            Lz = dims(3);  
            dz = Lz/(Nz-1);
            z = dz*(0:(Nz-1))' - Lz; % Cosine basis for DCT-I and DST-I
            nModes = Nz-1;
            
            if ~exist('rho0','var')
                rho0 = 1025;
            end
            rhoFunction = @(z) -(N0*N0*rho0/9.81)*z + rho0;
            rhobar = rhoFunction(z);
            N2 = N0*N0*ones(size(z));
            
            self@WaveVortexModel(dims, n, z, rhobar, N2, nModes, latitude, rho0);
            
            self.N0 = N0;
            
            self.BuildTransformationMatrices();
            
            % Preallocate this array for a faster dct
            self.realScratch = zeros(self.Nx,self.Ny,(2*self.Nz-1));
            self.complexScratch = complex(zeros(self.Nx,self.Ny,2*(self.Nz-1)));         
        end
                
        function h = get.h(self)
            [K,L,J] = ndgrid(self.k,self.l,self.j);
            K2 = K.*K + L.*L;
            M = J*pi/self.Lz;
            h = (1/self.g)*(self.N0*self.N0 - self.f0*self.f0)./(M.*M+K2);
            h(:,:,1) = 1; % prevent divide by zero
        end
                
        function self = BuildTransformationMatrices(self)  
            BuildTransformationMatrices@WaveVortexModel(self);
            
            % We renormalization the transformation matrices to directly
            % incorporate normalization of the modes and the DFT.          
            [~,~,J] = ndgrid(self.k,self.l,self.j);
            M = J*pi/self.Lz;
            N = self.N0;
            f = self.f0; 
            g_ = 9.81;
       
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Normalization for the vertical modes
            % This comes from equations B12 in the manuscript.
            signNorm = -2*(mod(J,2) == 1)+1; % equivalent to (-1)^j
            self.F = signNorm .* ((self.h).*M)*sqrt(2*g_/(self.Lz*(N*N-f*f)));
            self.G = signNorm .* sqrt(2*g_/(self.Lz*(N*N-f*f)));
            self.F(:,:,1) = 2; % j=0 mode is a factor of 2 too big in DCT-I
            self.G(:,:,1) = 1; % j=0 mode doesn't exist for G
  
            % Now make the Hermitian conjugate match.
            iFTransformScaling = 2./(self.Nx*self.Ny*self.F);
            iGTransformScaling = 2./(self.Nx*self.Ny*self.G);
            self.ApU = iFTransformScaling .* self.ApU;
            self.ApV = iFTransformScaling .* self.ApV;
            self.ApN = iGTransformScaling .* self.ApN;
            
            self.AmU = iFTransformScaling .* self.AmU;
            self.AmV = iFTransformScaling .* self.AmV;
            self.AmN = iGTransformScaling .* self.AmN;
            
            self.A0U = iFTransformScaling .* self.A0U;
            self.A0V = iFTransformScaling .* self.A0V;
            self.A0N = iGTransformScaling .* self.A0N;
                        
            % Now make the Hermitian conjugate match AND pre-multiply the
            % coefficients for the transformations.
            FTransformScaling = 0.5*self.Nx*self.Ny*self.F;
            self.UAp = FTransformScaling .* self.UAp;
            self.UAm = FTransformScaling .* self.UAm;
            self.UA0 = FTransformScaling .* self.UA0;
            
            self.VAp = FTransformScaling .* self.VAp;
            self.VAm = FTransformScaling .* self.VAm;
            self.VA0 = FTransformScaling .* self.VA0;
            
            GTransformScaling = 0.5*self.Nx*self.Ny*self.G;
            self.WAp = GTransformScaling .* self.WAp;
            self.WAm = GTransformScaling .* self.WAm;
            
            self.NAp = GTransformScaling .* self.NAp;
            self.NAm = GTransformScaling .* self.NAm;
            self.NA0 = GTransformScaling .* self.NA0;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Energetics
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function value = get.Apm_TE_factor(self)
            value = self.h; % factor of 2 larger than in the manuscript
            value(:,:,1) = self.Lz;
        end
        
        function value = get.A0_HKE_factor(self)
            [K,L,J] = ndgrid(self.k,self.l,self.j);
            K2 = K.*K + L.*L;
            M = J*pi/self.Lz;
            
            % This comes from equation (3.10) in the manuscript, but using
            % the relation from equation A2b
            % omega = sqrt(self.g*h.*K2 + self.f0*self.f0);
            % value = (self.g/(self.f0*self.f0)) * (omega.*omega - self.f0*self.f0) .* (self.N0*self.N0 - omega.*omega) / (2 * (self.N0*self.N0 - self.f0*self.f0) );
            value = (self.g^3/(self.f0*self.f0)) * K2.*self.h.*self.h.*M.*M / (2 * (self.N0*self.N0 - self.f0*self.f0) ); % factor of 2 larger than in the manuscript
            value(:,:,1) = (self.g^2/(self.f0*self.f0)) * K2(:,:,1) * self.Lz/2;
        end
        function value = get.A0_PE_factor(self)
            value = self.g*self.N0*self.N0/(self.N0*self.N0-self.f0*self.f0)/2; % factor of 2 larger than in the manuscript
        end
        function value = get.A0_TE_factor(self)
            value = self.A0_HKE_factor + self.A0_PE_factor;
        end
          
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Transformations to and from the spatial domain
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        function u_bar = TransformFromSpatialDomainWithF(self, u)
            % All coefficients are subsumbed into the transform
            % coefficients ApU,ApV,etc.
%             self.dctScratch = ifft(cat(3,u,u(:,:,(self.Nz-1):-1:2)),2*(self.Nz-1),3);
%             u_bar = real(self.dctScratch(:,:,1:(self.Nz-1))); % include barotropic mode, but leave off the Nyquist.
%             u_bar = fft(fft(u_bar,self.Nx,1),self.Ny,2);
%             
            u = ifft(cat(3,u,flip(u,3)),2*(self.Nz-1),3,'symmetric');
            u_bar = fft(fft(u(:,:,1:(self.Nz-1)),self.Nx,1),self.Ny,2);
        end
        
        function w_bar = TransformFromSpatialDomainWithG(self, w)
            % df = 1/(2*(Nz-1)*dz)
            % nyquist = (Nz-2)*df
            w = ifft(cat(3,w,-w(:,:,(self.Nz-1):-1:2)),2*(self.Nz-1),3);
            w_bar = imag(w(:,:,1:(self.Nz-1)));
            w_bar = fft(fft(w_bar,self.Nx,1),self.Ny,2);
        end
        
        function u = TransformToSpatialDomainWithF(self, u_bar)
            % All coefficients are subsumbed into the transform
            % coefficients UAp,UAm,etc.
            u = ifft(ifft(u_bar,self.Nx,1),self.Ny,2,'symmetric');    
            
            % Re-order to convert to a DCT-I via FFT
            u = cat(3,u,zeros(self.Nx,self.Ny),flip(u,3));
            u = ifft( u,2*(self.Nz-1),3,'symmetric');
            u = u(:,:,1:self.Nz)*(2*self.Nz-2); % We do not incorporate this coefficient into UAp, etc, so that the transforms remain inverses
        end  
                
        function w = TransformToSpatialDomainWithG(self, w_bar )
            % All coefficients are subsumbed into the transform
            % coefficients NAp,NAm,etc.
            w = ifft(ifft(w_bar,self.Nx,1),self.Ny,2,'symmetric');
            
            % Re-order to convert to a DST-I via FFT
            w = sqrt(-1)*cat(3,-w,zeros(self.Nx,self.Ny),flip(w,3));
            w = ifft( w,2*(self.Nz-1),3,'symmetric');
            w = w(:,:,1:self.Nz)*(2*self.Nz-2); % We do not incorporate this coefficient into UAp, etc, so that the transforms remain inverses
        end
        
        function [u,ux,uy,uz] = TransformToSpatialDomainWithFAllDerivatives(self, u_bar)
            % All coefficients are subsumbed into the transform
            % coefficients UAp,UAm,etc.
            u = ifft(ifft(u_bar,self.Nx,1),self.Ny,2,'symmetric');    
            
            % Re-order to convert to a DCT-I via FFT
            self.realScratch = cat(3,u,zeros(self.Nx,self.Ny),flip(u,3));
            u = ifft( self.realScratch,2*(self.Nz-1),3,'symmetric');
            u = u(:,:,1:self.Nz)*(2*self.Nz-2);
            
            ux = ifft( sqrt(-1)*self.k.*fft(u,self.Nx,1), self.Nx, 1,'symmetric');
            uy = ifft( sqrt(-1)*shiftdim(self.l,-1).*fft(u,self.Ny,2), self.Ny, 2,'symmetric');
            
            % To take the derivative, we multiply, then sine transform back
            m = reshape(-pi*self.j/self.Lz,1,1,[]);
            % That would be the normal multiplier, but we want to get it
            % ready for a DST, so.. -i*m,0,i*flip(m)
            m = cat(3,-sqrt(-1)*m,0,sqrt(-1)*flip(m,3));
            uz = ifft( m.*self.realScratch,2*(self.Nz-1),3,'symmetric');
            uz = uz(:,:,1:self.Nz)*(2*self.Nz-2);
        end  
        
        function [w,wx,wy,wz] = TransformToSpatialDomainWithGAllDerivatives(self, w_bar )
            % All coefficients are subsumbed into the transform
            % coefficients NAp,NAm,etc.
            w = ifft(ifft(w_bar,self.Nx,1),self.Ny,2,'symmetric');
            
            % Re-order to convert to a DST-I via FFT
            self.complexScratch = sqrt(-1)*cat(3,-w,zeros(self.Nx,self.Ny),flip(w,3));
            w = ifft( self.complexScratch,2*(self.Nz-1),3,'symmetric');
            w = w(:,:,1:self.Nz)*(2*self.Nz-2); % We do not incorporate this coefficient into UAp, etc, so that the transforms remain inverses
            
            wx = ifft( sqrt(-1)*self.k.*fft(w,self.Nx,1), self.Nx, 1,'symmetric');
            wy = ifft( sqrt(-1)*shiftdim(self.l,-1).*fft(w,self.Ny,2), self.Ny, 2,'symmetric');
            
            % To take the derivative, we multiply, then cosine transform back
            m = reshape(pi*self.j/self.Lz,1,1,[]);
            % That would be the normal multiplier, but we want to get it
            % ready for a DST, so.. i*m,0,-i*flip(m)
            m = cat(3,sqrt(-1)*m,0,-sqrt(-1)*flip(m,3));
            wz = ifft( m.*self.complexScratch,2*(self.Nz-1),3,'symmetric');
            wz = wz(:,:,1:self.Nz)*(2*self.Nz-2);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Needed to add and remove internal waves from the model
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function ratio = UmaxGNormRatioForWave(self,k0, l0, j0)
            if j0 == 0
                ratio = 1;
            else
                ratio = abs(1/self.F(k0+1,l0+1,j0+1));
            end
        end   
        
    end
   
        
        
end 



