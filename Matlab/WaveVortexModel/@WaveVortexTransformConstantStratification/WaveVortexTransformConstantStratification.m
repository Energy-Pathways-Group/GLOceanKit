classdef WaveVortexTransformConstantStratification < WaveVortexTransform
    %3D Boussinesq model with constant stratification solved in wave-vortex
    %space
    
    properties
        N0   
        realScratch, complexScratch; % of size Nx x Ny x (2*Nz-1)
        F,G
        h

        DCT, iDCT, DST, iDST, DFT, iDFT
        
        isHydrostatic = 0
        cg_x, cg_y, cg_z
        
        Apm_TE_factor
        A0_HKE_factor
        A0_PE_factor
        A0_TE_factor
    end
        
    methods
        function self = WaveVortexTransformConstantStratification(dims, n, latitude, N0, rho0, varargin)
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
            
            nargs = length(varargin);
            if mod(nargs,2) ~= 0
                error('Arguments must be given as name/value pairs');
            end
            isHydrostatic = 0;
            for k = 1:2:length(varargin)
                if strcmp(varargin{k}, 'hydrostatic')
                    isHydrostatic = varargin{k+1};
                end
            end

            Lz = dims(3);  
            dz = Lz/(Nz-1);
            z = dz*(0:(Nz-1))' - Lz; % Cosine basis for DCT-I and DST-I
            nModes = Nz-1;
            
            if ~exist('rho0','var') || isempty(rho0)
                rho0 = 1025;
            end
            rhoFunction = @(z) -(N0*N0*rho0/9.81)*z + rho0;
            N2Function = @(z) N0*N0*ones(size(z));
            rhobar = rhoFunction(z);
            N2 = N0*N0*ones(size(z));
            dLnN2 = zeros(size(z));
            
            self@WaveVortexTransform(dims, n, z, rhobar, N2, dLnN2, nModes, latitude, rho0);
            
            self.isHydrostatic = isHydrostatic;
            self.N0 = N0;
            self.rhoFunction = rhoFunction;
            self.N2Function = N2Function;
            
            self.buildTransformationMatrices();
            internalModes = InternalModesConstantStratification([N0 self.rho0], [-dims(3) 0],z,latitude);
            self.offgridModes = WaveVortexModelOffGrid(internalModes,latitude, N2Function,self.isHydrostatic);
            
            % Preallocate this array for a faster dct
            self.realScratch = zeros(self.Nx,self.Ny,(2*self.Nz-1));
            self.complexScratch = complex(zeros(self.Nx,self.Ny,2*(self.Nz-1)));
            warning('Need to check 2*(Nz-1), it gets extended to 2*Nz-1 during simulation');

            if 1 == 1
                self.DCT = CosineTransformForwardMatrix(self.Nz);
                self.DCT = self.DCT(1:end-1,:); % dump the Nyquist mode
                self.iDCT = CosineTransformBackMatrix(self.Nz);
                self.iDCT = self.iDCT(:,1:end-1); % dump the Nyquist mode
                self.DST = SineTransformForwardMatrix(self.Nz);
                self.DST = cat(1,zeros(1,self.Nz),self.DST);
                self.iDST = SineTransformBackMatrix(self.Nz);
                self.iDST = cat(2,zeros(self.Nz,1),self.iDST);
            end

        end
                
        function h = get.h(self)
            [K,L,J] = ndgrid(self.k,self.l,self.j);
            M = J*pi/self.Lz;
            if self.isHydrostatic == 1
                h = (1/self.g)*(self.N0*self.N0)./(M.*M);
                h(:,:,1) = 1; % prevent divide by zero
            else
                K2 = K.*K + L.*L;
                h = (1/self.g)*(self.N0*self.N0 - self.f0*self.f0)./(M.*M+K2);
                h(:,:,1) = 1; % prevent divide by zero
            end
        end
                
        function self = buildTransformationMatrices(self)

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
            if self.isHydrostatic == 1
                self.F = signNorm .* ((self.h).*M)*sqrt(2*g_/(self.Lz*N*N));
                self.G = signNorm .* sqrt(2*g_/(self.Lz*N*N));
            else
                self.F = signNorm .* ((self.h).*M)*sqrt(2*g_/(self.Lz*(N*N-f*f)));
                self.G = signNorm .* sqrt(2*g_/(self.Lz*(N*N-f*f)));
            end
            self.F(:,:,1) = 2; % j=0 mode is a factor of 2 too big in DCT-I
            self.G(:,:,1) = 1; % j=0 mode doesn't exist for G

            PP = self.F;
            QQ = self.G;

            buildTransformationMatrices@WaveVortexTransform(self,PP,QQ);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Nonlinear Flux
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [Fp,Fm,F0] = nonlinearFlux(self)
            uNL = self.u .* self.diffX(self.u)   + self.v .* self.diffY(self.u)   + self.w .* self.diffZF(self.u);
            vNL = self.u .* self.diffX(self.v)   + self.v .* self.diffY(self.v)   + self.w .* self.diffZF(self.v);
            nNL = self.u .* self.diffX(self.eta) + self.v .* self.diffY(self.eta) + self.w .* self.diffZG(self.eta);
            [Fp,Fm,F0] = wvt.transformUVEtaToWaveVortex(uNL,vNL,nNL,self.t);
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

            if self.isHydrostatic == 1
                value = (self.g^2/(self.f0*self.f0)) * K2 .* self.Apm_TE_factor/2;
            else
                M = J*pi/self.Lz;

                % This comes from equation (3.10) in the manuscript, but using
                % the relation from equation A2b
                % omega = sqrt(self.g*h.*K2 + self.f0*self.f0);
                % value = (self.g/(self.f0*self.f0)) * (omega.*omega - self.f0*self.f0) .* (self.N0*self.N0 - omega.*omega) / (2 * (self.N0*self.N0 - self.f0*self.f0) );
                value = (self.g^3/(self.f0*self.f0)) * K2.*self.h.*self.h.*M.*M / (2 * (self.N0*self.N0 - self.f0*self.f0) ); % factor of 2 larger than in the manuscript
                value(:,:,1) = (self.g^2/(self.f0*self.f0)) * K2(:,:,1) * self.Lz/2;
            end
        end
        function value = get.A0_PE_factor(self)
            if self.isHydrostatic == 1
                value = self.g*ones(self.Nx,self.Ny,self.nModes)/2;
            else
                value = self.g*self.N0*self.N0/(self.N0*self.N0-self.f0*self.f0)/2; % factor of 2 larger than in the manuscript
            end
        end
        function value = get.A0_TE_factor(self)
            value = self.A0_HKE_factor + self.A0_PE_factor;
        end
          
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Wave properties
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function cg_x = get.cg_x(self)
            [K,L,J] = ndgrid(self.k,self.l,self.j);
            K2 = K.*K + L.*L;
            M = J*pi/self.Lz;
            Omega = sqrt( (self.N0*self.N0*K2+self.f0*self.f0*M.*M)./(K2+M.*M) );
            cg_x = (K./Omega) .*M.*M .* (self.N0*self.N0-self.f0*self.f0)./(M.*M+K2).^2;
            cg_x(isnan(cg_x)) = 0;
        end
        
        function cg_y = get.cg_y(self)
            [K,L,J] = ndgrid(self.k,self.l,self.j);
            K2 = K.*K + L.*L;
            M = J*pi/self.Lz;
            Omega = sqrt( (self.N0*self.N0*K2+self.f0*self.f0*M.*M)./(K2+M.*M) );
            cg_y = (L./Omega) .* M.*M .* (self.N0*self.N0-self.f0*self.f0)./(M.*M+K2).^2;
            cg_y(isnan(cg_y)) = 0;
        end
        
        function cg_z = get.cg_z(self)
            [K,L,J] = ndgrid(self.k,self.l,self.j);
            K2 = K.*K + L.*L;
            M = J*pi/self.Lz;
            Omega = sqrt( (self.N0*self.N0*K2+self.f0*self.f0*M.*M)./(K2+M.*M) );
            cg_z = -(M./Omega) .* K2 .* (self.N0*self.N0-self.f0*self.f0)./(M.*M+K2).^2;
            cg_z(isnan(cg_z)) = 0;
        end
        
        function ProfileTransforms(self)
            Ubar = self.UAp.*self.Ap + self.UAm.*self.Am + self.UA0.*self.A0;
            Nbar = self.NAp.*self.Ap + self.NAm.*self.Am + self.NA0.*self.A0;
            MM = zeros(6,1);
            FF = zeros(6,1);
            names = {'    Finv', '    Ginv', '       F', '       G', 'Finv-all', 'Ginv-all'};
            N = 3;
            
            tic;
            for i=1:N
                U = self.transformToSpatialDomainWithF_MM(Ubar);
            end
            MM(1) = toc;
            tic;
            for i=1:N
                ETA = self.transformToSpatialDomainWithG_MM(Nbar);
            end
            MM(2) = toc;
            tic;
            for i=1:N
                self.transformFromSpatialDomainWithF_MM(U);
            end
            MM(3) = toc;
            tic;
            for i=1:N
                self.transformFromSpatialDomainWithG_MM(ETA);
            end
            MM(4) = toc;
            tic;
            for i=1:N
                self.transformToSpatialDomainWithFAllDerivatives_MM(Ubar);
            end
            MM(5) = toc;
            tic;
            for i=1:N
                self.transformToSpatialDomainWithGAllDerivatives_MM(Nbar);
            end
            MM(6) = toc;

            tic;
            for i=1:N
                U = self.transformToSpatialDomainWithF_FFT(Ubar);
            end
            FF(1) = toc;
            tic;
            for i=1:N
                ETA = self.transformToSpatialDomainWithG_FFT(Nbar);
            end
            FF(2) = toc;
            tic;
            for i=1:N
                self.transformFromSpatialDomainWithF_FFT(U);
            end
            FF(3) = toc;
            tic;
            for i=1:N
                self.transformFromSpatialDomainWithG_FFT(ETA);
            end
            FF(4) = toc;
            tic;
            for i=1:N
                self.transformToSpatialDomainWithFAllDerivatives_FFT(Ubar);
            end
            FF(5) = toc;
            tic;
            for i=1:N
                self.transformToSpatialDomainWithGAllDerivatives_FFT(Nbar);
            end
            FF(6) = toc;
            
            fprintf('--------|-- MM (s) --|-- FFT (s) --|\n')
            for i=1:6
                fprintf('%s|   %.4f   |   %.4f    | ', names{i}, MM(i), FF(i))
                if (MM(i)<FF(i))
                    fprintf('MM is %.2f faster\n',FF(i)/MM(i));
                else
                    fprintf('FFT is %.2f faster\n',MM(i)/FF(i));
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Transformations to and from the spatial domain, using FFTs
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        function u_bar = transformFromSpatialDomainWithF(self, u)     
            u_bar = self.transformFromSpatialDomainWithF_MM(u);
        end
        
        function w_bar = transformFromSpatialDomainWithG(self, w)
            w_bar = self.transformFromSpatialDomainWithG_MM(w);
        end
        
        function u = transformToSpatialDomainWithF(self, u_bar)
            u = self.transformToSpatialDomainWithF_MM(u_bar);
        end  
                
        function w = transformToSpatialDomainWithG(self, w_bar )
            w = self.transformToSpatialDomainWithG_MM(w_bar );
        end
        
        function [u,ux,uy,uz] = transformToSpatialDomainWithFAllDerivatives(self, u_bar)
            [u,ux,uy,uz] = self.transformToSpatialDomainWithFAllDerivatives_MM(u_bar);
        end  
        
        function [w,wx,wy,wz] = transformToSpatialDomainWithGAllDerivatives(self, w_bar )
            [w,wx,wy,wz] = self.transformToSpatialDomainWithGAllDerivatives_MM(w_bar );
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Transformations to and from the spatial domain, using matrix
        % multiplication (MM)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        function u_bar = transformFromSpatialDomainWithF_MM(self, u)
            u = permute(u,[3 1 2]); % keep adjacent in memory
            u = reshape(u,self.Nz,[]);
            u_bar = self.DCT*u;
            u_bar = reshape(u_bar,self.nModes,self.Nx,self.Ny);
            u_bar = permute(u_bar,[2 3 1]);

            u_bar = fft(fft(u_bar,self.Nx,1),self.Ny,2)/(self.Nx*self.Ny);
        end
        
        function w_bar = transformFromSpatialDomainWithG_MM(self, w)
            % df = 1/(2*(Nz-1)*dz)
            % nyquist = (Nz-2)*df
            w = permute(w,[3 1 2]); % keep adjacent in memory
            w = reshape(w,self.Nz,[]);
            w_bar = self.DST*w;
            w_bar = reshape(w_bar,self.nModes,self.Nx,self.Ny);
            w_bar = permute(w_bar,[2 3 1]);
            w_bar = fft(fft(w_bar,self.Nx,1),self.Ny,2)/(self.Nx*self.Ny);
        end
        
        function u = transformToSpatialDomainWithF_MM(self, u_bar)
            % All coefficients are subsumbed into the transform
            % coefficients UAp,UAm,etc.
            u_bar = ifft(ifft(u_bar,self.Nx,1),self.Ny,2,'symmetric')*self.Nx*self.Ny;    
            u_bar = permute(u_bar,[3 1 2]); % keep adjacent in memory
            u_bar = reshape(u_bar,self.nModes,[]);
            u = self.iDCT*u_bar;
            u = reshape(u,self.Nz,self.Nx,self.Ny);
            u = permute(u,[2 3 1]);
        end  
                
        function w = transformToSpatialDomainWithG_MM(self, w_bar )
            % All coefficients are subsumbed into the transform
            % coefficients NAp,NAm,etc.
            w_bar = ifft(ifft(w_bar,self.Nx,1),self.Ny,2,'symmetric')*self.Nx*self.Ny;
            w_bar = permute(w_bar,[3 1 2]); % keep adjacent in memory
            w_bar = reshape(w_bar,self.nModes,[]);
            w = self.iDST*w_bar;
            w = reshape(w,self.Nz,self.Nx,self.Ny);
            w = permute(w,[2 3 1]);        
        end
        
        function [u,ux,uy,uz] = transformToSpatialDomainWithFAllDerivatives_MM(self, u_bar)
            % All coefficients are subsumbed into the transform
            % coefficients UAp,UAm,etc.
            u_bar = ifft(ifft(u_bar,self.Nx,1),self.Ny,2,'symmetric')*self.Nx*self.Ny;    
            
            u_bar = permute(u_bar,[3 1 2]); % keep adjacent in memory
            u_bar = reshape(u_bar,self.nModes,[]);
            u = self.iDCT*u_bar;
            u = reshape(u,self.Nz,self.Nx,self.Ny);
            u = permute(u,[2 3 1]);
            
            ux = ifft( sqrt(-1)*self.k.*fft(u,self.Nx,1), self.Nx, 1,'symmetric');
            uy = ifft( sqrt(-1)*shiftdim(self.l,-1).*fft(u,self.Ny,2), self.Ny, 2,'symmetric');
            
            % To take the derivative, we multiply, then sine transform back
            m = reshape(-pi*self.j/self.Lz,[],1);
            uz = self.iDST*(m.*u_bar);
            uz = reshape(uz,self.Nz,self.Nx,self.Ny);
            uz = permute(uz,[2 3 1]);
        end  
        
        function [w,wx,wy,wz] = transformToSpatialDomainWithGAllDerivatives_MM(self, w_bar )
            % All coefficients are subsumbed into the transform
            % coefficients NAp,NAm,etc.
            w_bar = ifft(ifft(w_bar,self.Nx,1),self.Ny,2,'symmetric')*self.Nx*self.Ny;
            
            w_bar = permute(w_bar,[3 1 2]); % keep adjacent in memory
            w_bar = reshape(w_bar,self.nModes,[]);
            w = self.iDST*w_bar;
            w = reshape(w,self.Nz,self.Nx,self.Ny);
            w = permute(w,[2 3 1]);    
            
            wx = ifft( sqrt(-1)*self.k.*fft(w,self.Nx,1), self.Nx, 1,'symmetric');
            wy = ifft( sqrt(-1)*shiftdim(self.l,-1).*fft(w,self.Ny,2), self.Ny, 2,'symmetric');
            
            % To take the derivative, we multiply, then cosine transform back
            m = reshape(pi*self.j/self.Lz,[],1);
            wz = self.iDCT*(m.*w_bar);
            wz = reshape(wz,self.Nz,self.Nx,self.Ny);
            wz = permute(wz,[2 3 1]);
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Transformations to and from the spatial domain, using FFTs
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        function u_bar = transformFromSpatialDomainWithF_FFT(self, u)     
            u = ifft(cat(3,u,flip(u,3)),2*(self.Nz-1),3,'symmetric');
            u_bar = fft(fft(u(:,:,1:(self.Nz-1)),self.Nx,1),self.Ny,2)/(0.5*self.Nx*self.Ny);
        end
        
        function w_bar = transformFromSpatialDomainWithG_FFT(self, w)
            % df = 1/(2*(Nz-1)*dz)
            % nyquist = (Nz-2)*df
            w = ifft(cat(3,w,-w(:,:,(self.Nz-1):-1:2)),2*(self.Nz-1),3);
            w_bar = imag(w(:,:,1:(self.Nz-1)));
            w_bar = fft(fft(w_bar,self.Nx,1),self.Ny,2)/(0.5*self.Nx*self.Ny);
        end
        
        function u = transformToSpatialDomainWithF_FFT(self, u_bar)
            % All coefficients are subsumbed into the transform
            % coefficients UAp,UAm,etc.
            u = ifft(ifft(u_bar,self.Nx,1),self.Ny,2,'symmetric');    
            
            % Re-order to convert to a DCT-I via FFT
            u = cat(3,u,zeros(self.Nx,self.Ny),flip(u,3));
            u = ifft( u,2*(self.Nz-1),3,'symmetric');
            u = u(:,:,1:self.Nz)*(2*self.Nz-2)*0.5*self.Nx*self.Ny; % We do not incorporate this coefficient into UAp, etc, so that the transforms remain inverses
        end  
                
        function w = transformToSpatialDomainWithG_FFT(self, w_bar )
            % All coefficients are subsumbed into the transform
            % coefficients NAp,NAm,etc.
            w = ifft(ifft(w_bar,self.Nx,1),self.Ny,2,'symmetric');
            
            % Re-order to convert to a DST-I via FFT
            w = sqrt(-1)*cat(3,-w,zeros(self.Nx,self.Ny),flip(w,3));
            w = ifft( w,2*(self.Nz-1),3,'symmetric');
            w = w(:,:,1:self.Nz)*(2*self.Nz-2)*0.5*self.Nx*self.Ny; % We do not incorporate this coefficient into UAp, etc, so that the transforms remain inverses
        end
        
        function [u,ux,uy,uz] = transformToSpatialDomainWithFAllDerivatives_FFT(self, u_bar)
            % All coefficients are subsumbed into the transform
            % coefficients UAp,UAm,etc.
            u = ifft(ifft(u_bar,self.Nx,1),self.Ny,2,'symmetric')*0.5*self.Nx*self.Ny;    
            
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
        
        function [w,wx,wy,wz] = transformToSpatialDomainWithGAllDerivatives_FFT(self, w_bar )
            % All coefficients are subsumbed into the transform
            % coefficients NAp,NAm,etc.
            w = ifft(ifft(w_bar,self.Nx,1),self.Ny,2,'symmetric')*0.5*self.Nx*self.Ny;
            
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
        
        function ratio = uMaxGNormRatioForWave(self,k0, l0, j0)
            if j0 == 0
                ratio = 1;
            else
                ratio = abs(1/self.F(k0+1,l0+1,j0+1));
            end
        end   
        
    end
   
        
        
end 



