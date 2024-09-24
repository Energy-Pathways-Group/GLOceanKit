classdef WVTransformConstantStratification < WVTransform & WVStratifiedFlow & WVInertialOscillationMethods & WVGeostrophicMethods & WVMeanDensityAnomalyMethods & WVInternalGravityWaveMethods
    % Wave-vortex transformation that assumes constant stratification
    %
    % To initialization an instance of the
    % WVTransformConstantStratification class you must specific the
    % stratification and latitude, otherwise defaults will be assumed for
    % you.
    %
    % ```matlab
    % N0 = 3*2*pi/3600;
    % wvt = WVTransformConstantStratification([100e3, 100e3, 1300],[64, 64, 65], N0=N0,latitude=30);
    % ```
    %
    % - Topic: Initialization
    %
    % - Declaration: classdef WVTransformConstantStratification < [WVTransform](/classes/wvtransform/)
    properties (Access=protected) %(GetAccess=public, SetAccess=protected)
        F_g,G_g
        F_wg, G_wg

        DCT, iDCT, DST, iDST, DFT, iDFT

        cg_x, cg_y, cg_z
    end

    properties (GetAccess=public)
        N0
        h_pm
        h_0
        iOmega
        isHydrostatic = 0
    end
    properties (Dependent)
        FinvMatrix
        GinvMatrix
        FMatrix
        GMatrix
    end

    properties (Access=private)
        realScratch, complexScratch; % of size Nx x Ny x (2*Nz-1)
    end

    methods
        function self = WVTransformConstantStratification(Lxyz, Nxyz, options)
            % initialze a wave-vortex transform with constant stratification
            %
            %
            % To initialization an instance of the
            % WVTransformConstantStratification class you must specific the
            % stratification and latitude, otherwise defaults will be
            % assumed for you.
            %
            % To initialize a domain of size (100 km, 100 km, 1300 m) with
            % a buoyancy frequency of 3 cycles per hour at latitude 30,
            % call
            %
            % ```matlab
            % N0 = 3*2*pi/3600;
            % wvt = WVTransformConstantStratification([100e3, 100e3, 1300],[64, 64, 65], N0=N0,latitude=30);
            % ```
            %
            %
            %
            % - Topic: Initialization
            % - Declaration: wvt = WVTransformConstantStratification(Lxyz, Nxyz, options)
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
                options.latitude (1,1) double {mustBeGreaterThanOrEqual(options.latitude,5),mustBeLessThanOrEqual(options.latitude,85)} = 33
                options.rho0 (1,1) double {mustBePositive} = 1025
                options.isHydrostatic double {mustBeMember(options.isHydrostatic,[0 1])} = 0
            end
            Nz = Nxyz(3);
            % if mod(log2(Nxyz(3)-1),1) == 0
            %     Nz = Nxyz(3);
            % else
            %     error('The vertical dimension must have 2^n or (2^n)+1 points. This is an artificial restriction.');
            % end

            Lz = Lxyz(3);
            dz = Lz/(Nz-1);
            z = dz*(0:(Nz-1))' - Lz; % Cosine basis for DCT-I and DST-I
            nModes = Nz-1;
            N0 = options.N0;
            rho0 = options.rho0;
            rho = @(z) -(N0*N0*rho0/9.81)*z + rho0;
            N2 = @(z) N0*N0*ones(size(z));
            dLnN2 = @(z) zeros(size(z));
            verticalModes = InternalModesConstantStratification(N0=N0, rho0=rho0, zIn=[-Lz 0], zOut=z, latitude=options.latitude);

            self@WVStratifiedFlow(Lz,z,rho=rho,N2=N2,dLnN2=dLnN2,latitude=options.latitude,verticalModes=verticalModes)

            self@WVTransform(Lxyz, Nxyz(1:2), z, latitude=options.latitude,rho0=options.rho0,Nj=nModes);
            self.isHydrostatic = options.isHydrostatic;
            self.N0 = N0;

            self.buildVerticalModeProjectionOperators();

            self.initializeStratifiedFlow();
            self.initializeGeostrophicComponent();
            self.initializeMeanDensityAnomalyComponent();
            self.initializeInternalGravityWaveComponent();
            self.initializeInertialOscillationComponent();

            self.addPropertyAnnotations(WVPropertyAnnotation('N0',{},'rad s^{-1}', 'buoyancy frequency of the no-motion density'));

            % self.offgridModes = WVOffGridTransform(verticalModes,self.latitude, @(z) N0*N0*ones(size(z)),self.isHydrostatic);

            % Preallocate this array for a faster dct
            self.realScratch = zeros(self.Nx,self.Ny,(2*self.Nz-1));
            self.complexScratch = complex(zeros(self.Nx,self.Ny,2*(self.Nz-1)));
            % warning('Need to check 2*(Nz-1), it gets extended to 2*Nz-1 during simulation');

            self.DCT = WVTransformConstantStratification.CosineTransformForwardMatrix(self.Nz);
            self.DCT = self.DCT(1:self.Nj,:); % dump the Nyquist mode
            self.iDCT = WVTransformConstantStratification.CosineTransformBackMatrix(self.Nz);
            self.iDCT = self.iDCT(:,1:self.Nj); % dump the Nyquist mode
            self.DST = WVTransformConstantStratification.SineTransformForwardMatrix(self.Nz);
            self.DST = cat(1,zeros(1,self.Nz),self.DST);
            self.DST = self.DST(1:self.Nj,:);
            self.iDST = WVTransformConstantStratification.SineTransformBackMatrix(self.Nz);
            self.iDST = cat(2,zeros(self.Nz,1),self.iDST);
            self.iDST = self.iDST(:,1:self.Nj);

            %
            % outputVar = WVVariableAnnotation('rho_prime',{'x','y','z'},'kg/m3', 'density anomaly');
            % f = @(wvt) (wvt.rho0/9.81)*reshape(wvt.N2,1,1,[]).*wvt.transformToSpatialDomainWithG(wvt.NAp.*wvt.Apt + self.NAm.*wvt.Amt + self.NA0.*wvt.A0t);
            % self.addOperation(WVOperation('rho_prime',outputVar,f));
            %
            % outputVar = WVVariableAnnotation('rho_total',{'x','y','z'},'kg/m3', 'total potential density');
            % f = @(wvt) reshape(wvt.rhobar,1,1,[]) + wvt.rho_prime;
            % self.addOperation(WVOperation('rho_total',outputVar,f));

            self.nonlinearFluxOperation = WVNonlinearFlux(self);
        end

        function wvtX2 = waveVortexTransformWithResolution(self,m)
            wvtX2 = WVTransformConstantStratification([self.Lx self.Ly self.Lz],m, self.N0,latitude=self.latitude,rho0=self.rho0);
            wvtX2.t0 = self.t0;
            [wvtX2.Ap,wvtX2.Am,wvtX2.A0] = self.spectralVariableWithResolution(wvtX2,self.Ap,self.Am,self.A0);
            wvtX2.nonlinearFluxOperation = self.nonlinearFluxOperation.nonlinearFluxWithResolutionOfTransform(wvtX2);
        end

        function h = get.h_pm(self)
            [K,L,J] = self.kljGrid;
            M = J*pi/self.Lz;
            if self.isHydrostatic == 1
                h = (1/self.g)*(self.N0*self.N0)./(M.*M);
                h(J==0) = 1; % prevent divide by zero
            else
                K2 = K.*K + L.*L;
                h = (1/self.g)*(self.N0*self.N0 - self.f*self.f)./(M.*M+K2);
                h(J==0) = 1; % prevent divide by zero
            end
        end

        function h = get.h_0(self)
            M = reshape(self.j,[],1)*pi/self.Lz;
            h = (1/self.g)*(self.N0*self.N0)./(M.*M);
            h(1) = self.Lz;
        end

        function self = buildVerticalModeProjectionOperators(self)
            % We renormalization the transformation matrices to directly
            % incorporate normalization of the modes and the DFT.
            [~,~,J] = self.kljGrid;
            M = J*pi/self.Lz;
            N = self.N0;
            f = self.f;
            g_ = 9.81;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Normalization for the vertical modes
            % This comes from equations B12 in the manuscript.
            signNorm = -2*(mod(J,2) == 1)+1; % equivalent to (-1)^j
            self.F_g = signNorm .* ((self.h_0).*M)*sqrt(2*g_/(self.Lz*N*N));
            self.G_g = signNorm .* sqrt(2*g_/(self.Lz*N*N));
            self.F_g(J==0) = 2; % j=0 mode is a factor of 2 too big in DCT-I
            self.G_g(J==0) = 1; % j=0 mode doesn't exist for G

            if self.isHydrostatic == 1
                F_w = self.F_g;
                G_w = self.G_g;
            else
                F_w = signNorm .* ((self.h_pm).*M) * sqrt(2*g_/(self.Lz*(N*N-f*f)));
                G_w = signNorm .* sqrt(2*g_/(self.Lz*(N*N-f*f)));
            end
            F_w(J==0) = 2; % j=0 mode is a factor of 2 too big in DCT-I
            G_w(J==0) = 1;

            self.G_wg = self.G_g ./ G_w;
            self.F_wg = self.F_g ./ F_w;
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

        function u = transformToSpatialDomainWithF(self, options)
            arguments
                self WVTransform {mustBeNonempty}
                options.Apm double = 0
                options.A0 double = 0
            end
            if isscalar(options.Apm) && isscalar(options.A0)
                u = zeros(self.spatialMatrixSize);
            else
                u_bar = self.transformToSpatialDomainWithF_MM(options.Apm./self.F_wg + options.A0);
                u = self.transformToSpatialDomainWithFourier(u_bar);
            end
        end

        function w = transformToSpatialDomainWithG(self, options )
            arguments
                self WVTransform {mustBeNonempty}
                options.Apm double = 0
                options.A0 double = 0
            end
            if isscalar(options.Apm) && isscalar(options.A0)
                w = zeros(self.spatialMatrixSize);
            else
                w_bar = self.transformToSpatialDomainWithG_MM(options.Apm./self.G_wg + options.A0 );
                w = self.transformToSpatialDomainWithFourier(w_bar);
            end
        end


        % function [u,ux,uy,uz] = transformToSpatialDomainWithFAllDerivatives(self, options)
        %     arguments
        %         self WVTransform {mustBeNonempty}
        %         options.Apm double = 0
        %         options.A0 double = 0
        %     end
        %     [u,ux,uy,uz] = self.transformToSpatialDomainWithFAllDerivatives_MM(options.Apm./self.F_wg + options.A0);
        % end
        %
        % function [w,wx,wy,wz] = transformToSpatialDomainWithGAllDerivatives(self, options )
        %     arguments
        %         self WVTransform {mustBeNonempty}
        %         options.Apm double = 0
        %         options.A0 double = 0
        %     end
        %     [w,wx,wy,wz] = self.transformToSpatialDomainWithGAllDerivatives_MM(options.Apm./self.G_wg + options.A0 );
        % end


        function w_bar = transformWithG_wg(self, w_bar )
            w_bar = self.G_wg .* w_bar;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Transformations to and from the spatial domain, using matrix
        % multiplication (MM)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function u_bar = transformFromSpatialDomainWithFio(self,u)
            u_bar = self.DCT*u;
            u_bar = self.F_wg(:,1).*(u_bar./self.F_g(:,1));
        end

        function u_bar = transformFromSpatialDomainWithFg(self,u)
            u_bar = (self.DCT*u)./self.F_g;
        end

        function w_bar = transformFromSpatialDomainWithGg(self,w)
            w_bar = (self.DST*w)./self.G_g;
        end



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Transformation matrices
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function Finv = get.FinvMatrix(wvt)
            % transformation matrix $$F^{-1}$$
            %
            % A matrix that transforms a vector from vertical mode space to physical
            % space.
            %
            % - Topic: Operations — Transformations
            % - Declaration: Finv = FinvMatrix(wvt)
            % - Returns Finv: A matrix with dimensions [Nz Nj]
            arguments
                wvt         WVTransform
            end

            Finv = shiftdim(wvt.F_g(:,1),1) .* wvt.iDCT;

        end

        function F = get.FMatrix(wvt)
            % transformation matrix $$F$$
            %
            % A matrix that transforms a vector from physical
            % space to vertical mode space.
            %
            % - Topic: Operations — Transformations
            % - Declaration: F = FMatrix(wvt)
            % - Returns Finv: A matrix with dimensions [Nz Nj]
            arguments
                wvt         WVTransform
            end

            F = wvt.DCT ./ wvt.F_g(:,1);

        end
        function Ginv = get.GinvMatrix(wvt)
            % transformation matrix $$G^{-1}$$
            %
            % A matrix that transforms a vector from vertical mode space to physical
            % space.
            %
            % - Topic: Operations — Transformations
            % - Declaration: Ginv = GinvMatrix(wvt)
            % - Returns Finv: A matrix with dimensions [Nz Nj]
            arguments
                wvt         WVTransform
            end

            Ginv = shiftdim(wvt.G_g(:,1),1) .* wvt.iDST;

        end
        function G = get.GMatrix(wvt)
            % transformation matrix $$G$$
            %
            % A matrix that transforms a vector from physical
            % space to vertical mode space.
            %
            % - Topic: Operations — Transformations
            % - Declaration: G = GMatrix(wvt)
            % - Returns Ginv: A matrix with dimensions [Nz Nj]
            arguments
                wvt         WVTransform
            end

            G = wvt.DST ./ wvt.G_g(:,1);

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Needed to add and remove internal waves from the model
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function ratio = maxFw(self,kMode,lMode,j)
            arguments
                self WVTransform {mustBeNonempty}
                kMode (:,1) double
                lMode (:,1) double
                j (:,1) double
            end
            if j == 0
                ratio = 1;
            else
                indices = self.indexFromModeNumber(kMode,lMode,j);
                ratio = abs(self.F_g(indices)./self.F_wg(indices));
            end
        end
        function ratio = maxFg(self,kMode,lMode,j)
            arguments
                self WVTransform {mustBeNonempty}
                kMode (:,1) double
                lMode (:,1) double
                j (:,1) double
            end
            if j == 0
                ratio = 1;
            else
                indices = self.indexFromModeNumber(kMode,lMode,j);
                ratio = abs(self.F_g(indices));
            end
        end

        [ncfile,matFilePath] = writeToFile(wvt,path,variables,options)

        function flag = isequal(self,other)
            arguments
                self WVTransform
                other WVTransform
            end
            flag = isequal@WVTransform(self,other);
            flag = flag & isequal(self.N0, other.N0);
        end
    end



    methods (Static)
        wvt = waveVortexTransformFromFile(path,options)

        resultsTable = speedTest

        matrix = CosineTransformForwardMatrix(N)
        matrix = CosineTransformBackMatrix(n)
        matrix = SineTransformForwardMatrix(N)
        matrix = SineTransformBackMatrix(N)
    end


    methods %(Access=protected)
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

        function u_bar = transformFromSpatialDomainWithF_MM(self, u)
            u_bar = (self.DCT*u)./self.F_g;
        end

        function w_bar = transformFromSpatialDomainWithG_MM(self, w)
            w_bar = (self.DST*w)./self.G_g;
        end

        function u = transformToSpatialDomainWithF_MM(self, u_bar)
            u = self.iDCT*(self.F_g .* u_bar);
        end

        function w = transformToSpatialDomainWithG_MM(self, w_bar )
            w = self.iDST*(self.G_g .* w_bar);
        end

        function [u,ux,uy,uz] = transformToSpatialDomainWithFAllDerivatives_MM(self, u_bar)
            u_bar = self.F_g .* u_bar;

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
            w_bar = self.G_g .* w_bar;

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
            u_bar = (u_bar./self.F_g);
        end

        function w_bar = transformFromSpatialDomainWithG_FFT(self, w)
            % df = 1/(2*(Nz-1)*dz)
            % nyquist = (Nz-2)*df
            w = ifft(cat(3,w,-w(:,:,(self.Nz-1):-1:2)),2*(self.Nz-1),3);
            w_bar = imag(w(:,:,1:(self.Nz-1)));
            w_bar = fft(fft(w_bar,self.Nx,1),self.Ny,2)/(0.5*self.Nx*self.Ny);
            w_bar = (w_bar./self.G_g);
        end

        function u = transformToSpatialDomainWithF_FFT(self, u_bar)
            u_bar = self.F_g .* u_bar;

            % All coefficients are subsumbed into the transform
            % coefficients UAp,UAm,etc.
            u = ifft(ifft(u_bar,self.Nx,1),self.Ny,2,'symmetric');

            % Re-order to convert to a DCT-I via FFT
            u = cat(3,u,zeros(self.Nx,self.Ny),flip(u,3));
            u = ifft( u,2*(self.Nz-1),3,'symmetric');
            u = u(:,:,1:self.Nz)*(2*self.Nz-2)*0.5*self.Nx*self.Ny; % We do not incorporate this coefficient into UAp, etc, so that the transforms remain inverses
        end

        function w = transformToSpatialDomainWithG_FFT(self, w_bar )
            w_bar = self.G_g .* w_bar;

            % All coefficients are subsumbed into the transform
            % coefficients NAp,NAm,etc.
            w = ifft(ifft(w_bar,self.Nx,1),self.Ny,2,'symmetric');

            % Re-order to convert to a DST-I via FFT
            w = sqrt(-1)*cat(3,-w,zeros(self.Nx,self.Ny),flip(w,3));
            w = ifft( w,2*(self.Nz-1),3,'symmetric');
            w = w(:,:,1:self.Nz)*(2*self.Nz-2)*0.5*self.Nx*self.Ny; % We do not incorporate this coefficient into UAp, etc, so that the transforms remain inverses
        end

        function [u,ux,uy,uz] = transformToSpatialDomainWithFAllDerivatives_FFT(self, u_bar)
            u_bar = self.F_g .* u_bar;

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
            w_bar = self.G_g .* w_bar;

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
    end
end



