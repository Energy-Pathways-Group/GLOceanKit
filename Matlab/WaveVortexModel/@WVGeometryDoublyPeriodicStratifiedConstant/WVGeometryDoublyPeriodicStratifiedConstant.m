classdef WVGeometryDoublyPeriodicStratifiedConstant < WVGeometryDoublyPeriodic & WVStratification & WVGeometryCartesianXYZ
    properties (Access=public) %(GetAccess=public, SetAccess=protected)
        N0
        h_0
        h_pm

        F_g,G_g
        F_wg, G_wg
        DCT, iDCT, DST, iDST
    end

    properties
        isHydrostatic
    end

    properties (Dependent)
        FinvMatrix
        GinvMatrix
        FMatrix
        GMatrix
        Lr2
    end
    
    methods
        function self = WVGeometryDoublyPeriodicStratifiedConstant(Lxyz, Nxyz, geomOptions, stratOptions, options)
            % create geometry for 2D barotropic flow
            %
            % ```matlab
            % Lxy = 50e3;
            % Nxy = 256;
            % wvt = Cartesian2DBarotropic([Lxy, Lxy], [Nxy, Nxy]);
            % ```
            %
            % - Topic: Initialization
            % - Declaration: wvt = Cartesian2DBarotropic(Lxyz, Nxyz, options)
            % - Parameter Lxy: length of the domain (in meters) in the two coordinate directions, e.g. [Lx Ly]
            % - Parameter Nxy: number of grid points in the two coordinate directions, e.g. [Nx Ny]
            % - Parameter shouldAntialias: (optional) whether or not to de-alias for quadratic multiplications
            % - Returns wvt: a new Cartesian2DBarotropic instance
            arguments
                Lxyz (1,3) double {mustBePositive}
                Nxyz (1,3) double {mustBePositive}
                geomOptions.shouldAntialias (1,1) logical = true
                stratOptions.Nj (1,1) double {mustBePositive}
                stratOptions.rho0 (1,1) double {mustBePositive} = 1025
                stratOptions.planetaryRadius (1,1) double = 6.371e6
                stratOptions.rotationRate (1,1) double = 7.2921E-5
                stratOptions.latitude (1,1) double = 33
                stratOptions.g (1,1) double = 9.81
                options.isHydrostatic logical = false
                options.N0 (1,1) double {mustBePositive} = 5.2e-3
            end

            Lz = Lxyz(3);
            Nz = Nxyz(3);
            dz = Lz/(Nz-1);
            
            if ~isfield(stratOptions,"Nj")
                maxNj = Nxyz(3)-1;
                if geomOptions.shouldAntialias == true && maxNj > 3
                    stratOptions.Nj = floor(2*maxNj/3);
                else
                    stratOptions.Nj = maxNj;
                end
            end

            stratOptions.z = dz*(0:(Nz-1))' - Lz; % Cosine basis for DCT-I and DST-I
            stratOptions.j = (0:(stratOptions.Nj-1))';
            stratOptions.rhoFunction = @(z) -(options.N0*options.N0*stratOptions.rho0/stratOptions.g)*stratOptions.z + stratOptions.rho0;
            stratOptions.N2Function = @(z) options.N0*options.N0*ones(size(stratOptions.z));

            statOptionCell = namedargs2cell(stratOptions);
            self@WVStratification(Lz,Nz,statOptionCell{:});

            optionCell = namedargs2cell(geomOptions);
            self@WVGeometryDoublyPeriodic(Lxyz(1:2),Nxyz(1:2),optionCell{:},Nz=Nxyz(3),shouldExcludeNyquist=true,shouldExludeConjugates=true,conjugateDimension=2);

            self.isHydrostatic = options.isHydrostatic;
            self.N0 = options.N0;
            self.z_int = dz*ones(Nz,1);
            self.z_int(1) = dz/2; self.z_int(end) = dz/2; 
            self.buildVerticalModeProjectionOperators();
        end
        
        du = diffZF(self,u,options)
        dw = diffZG(self,w,options)

        function Lr2 = get.Lr2(self)
            Lr2 = self.g*self.h_0/(self.f*self.f);
        end

        function vm = verticalModes(self)
            vm = InternalModesConstantStratification(N0=self.N0, rho0=self.rho0, zIn=[self.z(1) 0], zOut=self.z, latitude=self.latitude);
        end

        function self = buildVerticalModeProjectionOperators(self)
            % Build the transformation matrices
            self.DCT = self.CosineTransformForwardMatrix(self.Nz);
            self.DCT = self.DCT(1:self.Nj,:); % dump the Nyquist mode
            self.iDCT = self.CosineTransformBackMatrix(self.Nz);
            self.iDCT = self.iDCT(:,1:self.Nj); % dump the Nyquist mode
            self.DST = self.SineTransformForwardMatrix(self.Nz);
            self.DST = cat(1,zeros(1,self.Nz),self.DST);
            self.DST = self.DST(1:self.Nj,:);
            self.iDST = self.SineTransformBackMatrix(self.Nz);
            self.iDST = cat(2,zeros(self.Nz,1),self.iDST);
            self.iDST = self.iDST(:,1:self.Nj);

            % We renormalization the transformation matrices to directly
            % incorporate normalization of the modes and the DFT.
            [~,~,J] = self.kljGrid;
            M = J*pi/self.Lz;
            N = self.N0;
            g_ = 9.81;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Normalization for the vertical modes
            % This comes from equations B12 in the manuscript.
            signNorm = -2*(mod(J,2) == 1)+1; % equivalent to (-1)^j
            self.F_g = signNorm .* ((self.h_0).*M)*sqrt(2*g_/(self.Lz*N*N));
            self.G_g = signNorm .* sqrt(2*g_/(self.Lz*N*N));
            self.F_g(J==0) = 2; % j=0 mode is a factor of 2 too big in DCT-I
            self.G_g(J==0) = 1; % j=0 mode doesn't exist for G

            if self.isHydrostatic == true
                F_w = self.F_g;
                G_w = self.G_g;
            else
                f = self.f;
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
        % Transformation matrices (geostrophic)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function h = get.h_0(self)
            M = reshape(self.j,[],1)*pi/self.Lz;
            h = (1/self.g)*(self.N0*self.N0)./(M.*M);
            h(1) = self.Lz;
        end

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
        % Transformation matrices (wave)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function h = get.h_pm(self)
            [K,L,J] = self.kljGrid;
            M = J*pi/self.Lz;
            if self.isHydrostatic == true
                h = (1/self.g)*(self.N0*self.N0)./(M.*M);
                h(J==0) = 1; % prevent divide by zero
            else
                K2 = K.*K + L.*L;
                h = (1/self.g)*(self.N0*self.N0 - self.f*self.f)./(M.*M+K2);
                h(J==0) = 1; % prevent divide by zero
            end
        end

        function Finv = FwInvMatrix(self,kMode,lMode)
            % transformation matrix $$F_w^{-1}$$
            %
            % A matrix that transforms a vector of igw amplitudes from
            % vertical mode space to physical space.
            %
            % - Topic: Operations — Transformations
            % - Declaration: Finv = FwInvMatrix(wvt,kMode,lMode)
            % - Returns Finv: A matrix with dimensions [Nz Nj]
            arguments
                self WVTransform
                kMode (1,1) double
                lMode (1,1) double
            end

            iK = self.indexFromModeNumber(kMode,lMode,self.j);
            Finv = shiftdim(self.F_g(iK)./self.F_wg(iK),1) .* self.iDCT; % 
        end

        function F = FwMatrix(self,kMode,lMode)
            % transformation matrix $$F_w$$
            %
            % A matrix that transforms a vector in physical space to IGW
            % mode space
            %
            % - Topic: Operations — Transformations
            % - Declaration: F = FwMatrix(wvt,kMode,lMode)
            % - Returns F: A matrix with dimensions [Nj Nz]
            arguments
                self WVTransform
                kMode (1,1) double
                lMode (1,1) double
            end

            iK = self.indexFromModeNumber(kMode,lMode,self.j);
            F = self.DCT ./ reshape(self.F_g(iK)./self.F_wg(iK),[],1);
        end

        function Ginv = GwInvMatrix(self,kMode,lMode)
            % transformation matrix $$G_w^{-1}$$
            %
            % A matrix that transforms a vector of igw amplitudes from
            % vertical mode space to physical space.
            %
            % - Topic: Operations — Transformations
            % - Declaration: Ginv = GwInvMatrix(wvt,kMode,lMode)
            % - Returns Ginv: A matrix with dimensions [Nz Nj]
            arguments
                self WVTransform
                kMode (1,1) double
                lMode (1,1) double
            end

            iK = self.indexFromModeNumber(kMode,lMode,self.j);
            Ginv = shiftdim(self.G_g(iK)./self.G_wg(iK),1) .* self.iDST; %
        end

        function G = GwMatrix(self,kMode,lMode)
            % transformation matrix $$G_w$$
            %
            % A matrix that transforms a vector in physical space to IGW
            % mode space
            %
            % - Topic: Operations — Transformations
            % - Declaration: G = GwMatrix(wvt,kMode,lMode)
            % - Returns G: A matrix with dimensions [Nj Nz]
            arguments
                self WVTransform
                kMode (1,1) double
                lMode (1,1) double
            end

            iK = self.indexFromModeNumber(kMode,lMode,self.j);
            G = self.DST ./ reshape(self.G_g(iK)./self.G_wg(iK),[],1);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Norm ratio (needed to add and remove internal waves from the
        % model)
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Transformation functions
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

        function w_bar = transformWithG_wg(self, w_bar )
            w_bar = self.G_wg .* w_bar;
        end

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
    end

    methods (Access=protected)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Transformations to and from the spatial domain, using matrix
        % multiplication (MM)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Profile the different transformations
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    end

    methods (Static)
        matrix = CosineTransformForwardMatrix(N)
        matrix = CosineTransformBackMatrix(n)
        matrix = SineTransformForwardMatrix(N)
        matrix = SineTransformBackMatrix(N)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % CAAnnotatedClass required methods, which enables writeToFile
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = WVGeometryDoublyPeriodicStratifiedConstant.propertyAnnotationsForGeometry();
        end

        function vars = classRequiredPropertyNames()
            vars = WVGeometryDoublyPeriodicStratifiedConstant.namesOfRequiredPropertiesForGeometry();
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Stratification specific property annotations and initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function requiredPropertyNames = namesOfRequiredPropertiesForGeometry()
            requiredPropertyNames = WVStratification.namesOfRequiredPropertiesForStratification();
            requiredPropertyNames = union(requiredPropertyNames,WVGeometryDoublyPeriodic.namesOfRequiredPropertiesForGeometry());
            requiredPropertyNames = union(requiredPropertyNames,WVGeometryDoublyPeriodicStratifiedConstant.newRequiredPropertyNames());
            requiredPropertyNames = setdiff(requiredPropertyNames,WVGeometryDoublyPeriodicStratifiedConstant.newNonrequiredPropertyNames);
        end

        function newNonrequiredPropertyNames = newNonrequiredPropertyNames()
            newNonrequiredPropertyNames = {'N2Function','j','conjugateDimension','shouldExcludeNyquist','shouldExludeConjugates'};
        end

        function newRequiredPropertyNames = newRequiredPropertyNames()
            newRequiredPropertyNames = {'isHydrostatic', 'N0'};
        end

        function propertyAnnotations = propertyAnnotationsForGeometry()
            propertyAnnotations = WVGeometryDoublyPeriodic.propertyAnnotationsForGeometry();
            propertyAnnotations = cat(2,propertyAnnotations,WVStratification.propertyAnnotationsForStratification());
            propertyAnnotations = cat(2,propertyAnnotations,WVGeometryCartesianXYZ.propertyAnnotationsForGeometry());
            propertyAnnotations(end+1) = CANumericProperty('isHydrostatic',{},'bool', 'whether the transforms are hydrostatic or non-hydrostatic', detailedDescription='- topic: Domain Attributes — Grid');
            propertyAnnotations(end+1) = CANumericProperty('N0',{},'rad s^{-1}', 'buoyancy frequency of the no-motion density', detailedDescription="The buoyancy frequency $$N_0$$ is defined as $$N_0\equiv sqrt{ - \frac{g}{\rho_0} \frac{\partial \rho_\textrm{nm}}{\partial z} }$$.");
            propertyAnnotations(end+1) = CANumericProperty('h_0',{'j'},'m', 'equivalent depth of each geostrophic mode', detailedDescription='- topic: Domain Attributes — Stratification');
            propertyAnnotations(end+1) = CANumericProperty('Lr2',{'j'},'m^2', 'squared Rossby radius');
        end

        function [Lxyz, Nxyz, options] = requiredPropertiesForGeometryFromGroup(group)
            arguments (Input)
                group NetCDFGroup {mustBeNonempty}
            end
            arguments (Output)
                Lxyz (1,3) double {mustBePositive}
                Nxyz (1,3) double {mustBePositive}
                options
            end
            [Lxyz(1:2), Nxyz(1:2), geomOptions] = WVGeometryDoublyPeriodic.requiredPropertiesForGeometryFromGroup(group,shouldIgnoreMissingProperties=true);
            [Lxyz(3), Nxyz(3), stratOptions] = WVStratification.requiredPropertiesForStratificationFromGroup(group,shouldIgnoreMissingProperties=true);
            vars = CAAnnotatedClass.propertyValuesFromGroup(group,WVGeometryDoublyPeriodicStratifiedConstant.newRequiredPropertyNames);
            S = struct(stratOptions{:});
            S = rmfield(S,'j');
            S = rmfield(S,'z');
            stratOptions = namedargs2cell(S);
            newOptions = namedargs2cell(vars);
            options = cat(2,stratOptions,geomOptions,newOptions);
        end

        function geometry = geometryFromFile(path)
            arguments (Input)
                path char {mustBeNonempty}
            end
            arguments (Output)
                geometry WVGeometryDoublyPeriodicStratifiedConstant {mustBeNonempty}
            end
            ncfile = NetCDFFile(path);
            geometry = WVGeometryDoublyPeriodicStratifiedConstant.geometryFromGroup(ncfile);
        end

        function geometry = geometryFromGroup(group)
            arguments (Input)
                group NetCDFGroup {mustBeNonempty}
            end
            arguments (Output)
                geometry WVGeometryDoublyPeriodicStratifiedConstant {mustBeNonempty}
            end
            CAAnnotatedClass.throwErrorIfMissingProperties(group,WVGeometryDoublyPeriodicStratifiedConstant.namesOfRequiredPropertiesForGeometry);
            [Lxyz, Nxyz, options] = WVGeometryDoublyPeriodicStratifiedConstant.requiredPropertiesForGeometryFromGroup(group);
            geometry = WVGeometryDoublyPeriodicStratifiedConstant(Lxyz,Nxyz,options{:});
        end

    end
end