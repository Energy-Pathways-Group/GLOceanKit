classdef WVTransformConstantStratification < WVTransform
    % Wave-vortex transformation that assumes constant stratification
    %
    % To initialization an instance of the
    % WVTransformConstantStratification class you must specific the
    % stratification and latitude, otherwise defaults will be assumed for
    % you.
    % 
    % ```matlab
    % N0 = 3*2*pi/3600;
    % wvt = WVTransformConstantStratification([100e3, 100e3, 1300],[64, 64, 65], NN0=N0,latitude=30);
    % ```
    %
    % - Topic: Initialization
    %
    % - Declaration: classdef WVTransformConstantStratification < [WVTransform](/classes/wvtransform/)
    properties (GetAccess=public, SetAccess=protected)
        N0, N2, rhobar
        F,G
        h

        DCT, iDCT, DST, iDST, DFT, iDFT
        
        isHydrostatic = 0
        cg_x, cg_y, cg_z
        
        Apm_TE_factor
        A0_HKE_factor
        A0_PE_factor
        A0_TE_factor
        A0_TZ_factor
        A0_QGPV_factor
    end

    properties (Access=private)
        realScratch, complexScratch; % of size Nx x Ny x (2*Nz-1)
    end
        
    methods
        function self = WVTransformConstantStratification(Lxyz, Nxyz, options)
            % initialze a wave-vortex transform with constant stratification
            %
            % - Topic: Initialization
            % - Declaration: wvt = WVTransformConstantStratification(Lxyz, Nxyz, N0, options)
            % - Parameter Lxyz: length of the domain (in meters) in the three coordinate directions, e.g. [Lx Ly Lz]
            % - Parameter Nxyz: number of grid points in the three coordinate directions, e.g. [Nx Ny Nz]
            % - Parameter N0:  (optional) buoyancy frequency (radians/s) default is 5.2e-3, or 3 cph)
            % - Parameter latitude: (optional) latitude of the domain (default is 33 degrees north)
            % - Parameter rho0: (optional) density at the surface z=0 (default is 1025 kg/m^3)
            % - Parameter isHydrostatic: (optional) flag indicating whether to use hydrostatic transformations (default 0)
            % - Returns wvt: a new WVTransformConstantStratification instance
            arguments
                Lxyz (1,3) double {mustBePositive}
                Nxyz (1,3) double {mustBePositive}
                options.N0 (1,1) double {mustBePositive} = 5.2e-3
                options.latitude (1,1) double = 33
                options.rho0 (1,1) double {mustBePositive} = 1025
                options.isHydrostatic double {mustBeMember(options.isHydrostatic,[0 1])} = 0 
            end
                        
            if mod(log2(Nxyz(3)-1),1) == 0
                Nz = Nxyz(3);
            else
                error('The vertical dimension must have 2^n or (2^n)+1 points. This is an artificial restriction.');
            end

            Lz = Lxyz(3);  
            dz = Lz/(Nz-1);
            z = dz*(0:(Nz-1))' - Lz; % Cosine basis for DCT-I and DST-I
            nModes = Nz-1;
            N0 = options.N0;
            
            self@WVTransform(Lxyz, Nxyz(1:2), z, latitude=options.latitude,rho0=options.rho0,Nj=nModes,Nmax=N0);
            
            self.isHydrostatic = options.isHydrostatic;
            self.N0 = N0;
            self.N2 = N0*N0;
            rhoFunction = @(z) -(N0*N0*self.rho0/9.81)*z + self.rho0;
            self.rhobar = rhoFunction(z);

            self.buildTransformationMatrices();
%             internalModes = InternalModesConstantStratification([N0 self.rho0], [-Lxyz(3) 0],z,self.latitude);
            internalModes = InternalModesConstantStratification(N0=N0, rho0=self.rho0, zIn=[-Lxyz(3) 0], zOut=z, latitude=self.latitude);
            self.offgridModes = WVOffGridTransform(internalModes,self.latitude, @(z) N0*N0*ones(size(z)),self.isHydrostatic);
            
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

            % 
            % outputVar = WVVariableAnnotation('rho_prime',{'x','y','z'},'kg/m3', 'density anomaly');
            % f = @(wvt) (wvt.rho0/9.81)*reshape(wvt.N2,1,1,[]).*wvt.transformToSpatialDomainWithG(wvt.NAp.*wvt.Apt + self.NAm.*wvt.Amt + self.NA0.*wvt.A0t);
            % self.addOperation(WVOperation('rho_prime',outputVar,f));
            % 
            % outputVar = WVVariableAnnotation('rho_total',{'x','y','z'},'kg/m3', 'total potential density');
            % f = @(wvt) reshape(wvt.rhobar,1,1,[]) + wvt.rho_prime;
            % self.addOperation(WVOperation('rho_total',outputVar,f));

            self.nonlinearFluxOperation = Boussinesq(self);
        end
                
        function wvtX2 = waveVortexTransformWithResolution(self,m)
            wvtX2 = WVTransformConstantStratification([self.Lx self.Ly self.Lz],m, self.N0,latitude=self.latitude,rho0=self.rho0);
            wvtX2.t0 = self.t0;
            if wvtX2.Nx>=self.Nx && wvtX2.Ny >= self.Ny && wvtX2.Nj >= self.Nj
                kIndices = cat(2,1:(self.Nk/2),(wvtX2.Nk-self.Nk/2 + 1):wvtX2.Nk);
                lIndices = cat(2,1:(self.Nl/2),(wvtX2.Nl-self.Nl/2 + 1):wvtX2.Nl);
                wvtX2.Ap(kIndices,lIndices,1:self.Nj) = self.Ap;
                wvtX2.Am(kIndices,lIndices,1:self.Nj) = self.Am;
                wvtX2.A0(kIndices,lIndices,1:self.Nj) = self.A0;
            else
                error('Reducing resolution not yet implemented. Go for it though, it should be easy.');
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
                h = (1/self.g)*(self.N0*self.N0 - self.f*self.f)./(M.*M+K2);
                h(:,:,1) = 1; % prevent divide by zero
            end
        end
                
        function self = buildTransformationMatrices(self)

            % We renormalization the transformation matrices to directly
            % incorporate normalization of the modes and the DFT.          
            [~,~,J] = ndgrid(self.k,self.l,self.j);
            M = J*pi/self.Lz;
            N = self.N0;
            f = self.f; 
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

            buildTransformationMatrices@WVTransform(self);
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
                value = (self.g^2/(self.f*self.f)) * K2 .* self.Apm_TE_factor/2;
            else
                M = J*pi/self.Lz;

                % This comes from equation (3.10) in the manuscript, but using
                % the relation from equation A2b
                % omega = sqrt(self.g*h.*K2 + self.f*self.f);
                % value = (self.g/(self.f*self.f)) * (omega.*omega - self.f*self.f) .* (self.N0*self.N0 - omega.*omega) / (2 * (self.N0*self.N0 - self.f*self.f) );
                value = (self.g^3/(self.f*self.f)) * K2.*self.h.*self.h.*M.*M / (2 * (self.N0*self.N0 - self.f*self.f) ); % factor of 2 larger than in the manuscript
                value(:,:,1) = (self.g^2/(self.f*self.f)) * K2(:,:,1) * self.Lz/2;
            end
        end
        function value = get.A0_PE_factor(self)
            if self.isHydrostatic == 1
                value = self.g*ones(self.Nk,self.Nl,self.Nj)/2;
                value(:,:,1) = 0;
            else
                value = self.g*self.N0*self.N0/(self.N0*self.N0-self.f*self.f)/2; % factor of 2 larger than in the manuscript
                value(:,:,1) = 0;
            end
        end
        function value = get.A0_TE_factor(self)
            value = self.A0_HKE_factor + self.A0_PE_factor;
        end

        function value = get.A0_TZ_factor(self)
            error('Not yet implemented for constant stratification');
            % Kh = self.Kh;
            % Lr2 = wvt.g*(wvt.h)/(wvt.f*wvt.f);
            % Lr2(1) = wvt.g*wvt.Lz/(wvt.f*wvt.f);
            % value = (self.g/2) * Lr2 .* ( (self.Kh).^2 + Lr2.^(-1) ).^2;
            % value(:,:,1) = (self.g/2) * Lr2(1) .* (Kh(:,:,1)).^4;
        end

        function value = get.A0_QGPV_factor(self)
            error('Not yet implemented for constant stratification');
            % Kh = self.Kh;
            % Lr2 = wvt.g*(wvt.h)/(wvt.f*wvt.f);
            % Lr2(1) = wvt.g*wvt.Lz/(wvt.f*wvt.f);
            % value = (self.g/2) * Lr2 .* ( (self.Kh).^2 + Lr2.^(-1) ).^2;
            % value(:,:,1) = (self.g/2) * Lr2(1) .* (Kh(:,:,1)).^4;
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
            Omega = sqrt( (self.N0*self.N0*K2+self.f*self.f*M.*M)./(K2+M.*M) );
            cg_x = (K./Omega) .*M.*M .* (self.N0*self.N0-self.f*self.f)./(M.*M+K2).^2;
            cg_x(isnan(cg_x)) = 0;
        end
        
        function cg_y = get.cg_y(self)
            [K,L,J] = ndgrid(self.k,self.l,self.j);
            K2 = K.*K + L.*L;
            M = J*pi/self.Lz;
            Omega = sqrt( (self.N0*self.N0*K2+self.f*self.f*M.*M)./(K2+M.*M) );
            cg_y = (L./Omega) .* M.*M .* (self.N0*self.N0-self.f*self.f)./(M.*M+K2).^2;
            cg_y(isnan(cg_y)) = 0;
        end
        
        function cg_z = get.cg_z(self)
            [K,L,J] = ndgrid(self.k,self.l,self.j);
            K2 = K.*K + L.*L;
            M = J*pi/self.Lz;
            Omega = sqrt( (self.N0*self.N0*K2+self.f*self.f*M.*M)./(K2+M.*M) );
            cg_z = -(M./Omega) .* K2 .* (self.N0*self.N0-self.f*self.f)./(M.*M+K2).^2;
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
            u_bar = reshape(u_bar,self.Nj,self.Nx,self.Ny);
            u_bar = permute(u_bar,[2 3 1]);
            u_bar = fft(fft(u_bar,self.Nx,1),self.Ny,2);
            u_bar = (u_bar./self.F)/(self.Nx*self.Ny);
        end
        
        function w_bar = transformFromSpatialDomainWithG_MM(self, w)
            % df = 1/(2*(Nz-1)*dz)
            % nyquist = (Nz-2)*df
            w = permute(w,[3 1 2]); % keep adjacent in memory
            w = reshape(w,self.Nz,[]);
            w_bar = self.DST*w;
            w_bar = reshape(w_bar,self.Nj,self.Nx,self.Ny);
            w_bar = permute(w_bar,[2 3 1]);
            w_bar = fft(fft(w_bar,self.Nx,1),self.Ny,2);
            w_bar = (w_bar./self.G)/(self.Nx*self.Ny);
        end
        
        function u = transformToSpatialDomainWithF_MM(self, u_bar)
            % All coefficients are subsumbed into the transform
            % coefficients UAp,UAm,etc.
            u_bar = self.F .* u_bar;
            u_bar = ifft(ifft(u_bar,self.Nx,1),self.Ny,2,'symmetric')*self.Nx*self.Ny;    
            u_bar = permute(u_bar,[3 1 2]); % keep adjacent in memory
            u_bar = reshape(u_bar,self.Nj,[]);
            u = self.iDCT*u_bar;
            u = reshape(u,self.Nz,self.Nx,self.Ny);
            u = permute(u,[2 3 1]);
        end  
                
        function w = transformToSpatialDomainWithG_MM(self, w_bar )
            % All coefficients are subsumbed into the transform
            % coefficients NAp,NAm,etc.
            w_bar = self.G .* w_bar;
            w_bar = ifft(ifft(w_bar,self.Nx,1),self.Ny,2,'symmetric')*self.Nx*self.Ny;
            w_bar = permute(w_bar,[3 1 2]); % keep adjacent in memory
            w_bar = reshape(w_bar,self.Nj,[]);
            w = self.iDST*w_bar;
            w = reshape(w,self.Nz,self.Nx,self.Ny);
            w = permute(w,[2 3 1]);        
        end
        
        function [u,ux,uy,uz] = transformToSpatialDomainWithFAllDerivatives_MM(self, u_bar)
            u_bar = self.F .* u_bar;

            % All coefficients are subsumbed into the transform
            % coefficients UAp,UAm,etc.
            u_bar = ifft(ifft(u_bar,self.Nx,1),self.Ny,2,'symmetric')*self.Nx*self.Ny;    
            
            u_bar = permute(u_bar,[3 1 2]); % keep adjacent in memory
            u_bar = reshape(u_bar,self.Nj,[]);
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
            w_bar = self.G .* w_bar;

            % All coefficients are subsumbed into the transform
            % coefficients NAp,NAm,etc.
            w_bar = ifft(ifft(w_bar,self.Nx,1),self.Ny,2,'symmetric')*self.Nx*self.Ny;
            
            w_bar = permute(w_bar,[3 1 2]); % keep adjacent in memory
            w_bar = reshape(w_bar,self.Nj,[]);
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
            u_bar = (u_bar./self.F);
        end
        
        function w_bar = transformFromSpatialDomainWithG_FFT(self, w)
            % df = 1/(2*(Nz-1)*dz)
            % nyquist = (Nz-2)*df
            w = ifft(cat(3,w,-w(:,:,(self.Nz-1):-1:2)),2*(self.Nz-1),3);
            w_bar = imag(w(:,:,1:(self.Nz-1)));
            w_bar = fft(fft(w_bar,self.Nx,1),self.Ny,2)/(0.5*self.Nx*self.Ny);
            w_bar = (w_bar./self.G);
        end
        
        function u = transformToSpatialDomainWithF_FFT(self, u_bar)
            u_bar = self.F .* u_bar;

            % All coefficients are subsumbed into the transform
            % coefficients UAp,UAm,etc.
            u = ifft(ifft(u_bar,self.Nx,1),self.Ny,2,'symmetric');    
            
            % Re-order to convert to a DCT-I via FFT
            u = cat(3,u,zeros(self.Nx,self.Ny),flip(u,3));
            u = ifft( u,2*(self.Nz-1),3,'symmetric');
            u = u(:,:,1:self.Nz)*(2*self.Nz-2)*0.5*self.Nx*self.Ny; % We do not incorporate this coefficient into UAp, etc, so that the transforms remain inverses
        end  
                
        function w = transformToSpatialDomainWithG_FFT(self, w_bar )
            w_bar = self.G .* w_bar;

            % All coefficients are subsumbed into the transform
            % coefficients NAp,NAm,etc.
            w = ifft(ifft(w_bar,self.Nx,1),self.Ny,2,'symmetric');
            
            % Re-order to convert to a DST-I via FFT
            w = sqrt(-1)*cat(3,-w,zeros(self.Nx,self.Ny),flip(w,3));
            w = ifft( w,2*(self.Nz-1),3,'symmetric');
            w = w(:,:,1:self.Nz)*(2*self.Nz-2)*0.5*self.Nx*self.Ny; % We do not incorporate this coefficient into UAp, etc, so that the transforms remain inverses
        end
        
        function [u,ux,uy,uz] = transformToSpatialDomainWithFAllDerivatives_FFT(self, u_bar)
            u_bar = self.F .* u_bar;

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
            w_bar = self.G .* w_bar;

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
        
        [ncfile,matFilePath] = writeToFile(wvt,path,variables,options)
    end

    methods (Static)
        wvt = waveVortexTransformFromFile(path,options)
    end
        
end 



