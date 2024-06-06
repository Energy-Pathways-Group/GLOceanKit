classdef WVTransformBoussinesq < WVTransform & WVStratifiedFlow & WVInertialOscillationMethods & WVGeostrophicMethods & WVMeanDensityAnomalyMethods & WVInternalGravityWaveMethods
    % 3D hydrostatic Boussinesq model with arbitrary stratification solved
    % in wave-vortex space
    %
    % Couple of different initialization paths:
    % 1) You want to run this as a prognostic model and therefore want
    %    the chebyshev points automatically found for you
    %       Init([Lx Ly Lz], [Nx Ny Nz], latitude, rho)
    %
    % 2) You want to run this as a diagnostic model and therefore want
    %    to specify the depths and modes yourself
    %       Init([Lx Ly Lz], [Nx Ny Nz], latitude, rho, 'zgrid', z)

    properties (Access=protected) %(GetAccess=public, SetAccess=protected)
        % Geostrophic transformation matrices
        PF0inv, QG0inv % size(PFinv,PGinv)=[Nz x Nj x 1]
        PF0, QG0 % size(PF,PG)=[Nj x Nz x 1]
        P0 % Preconditioner for F, size(P)=[1 1 Nj x 1]. F*u = uhat, (PF)*u = P*uhat, so ubar==P*uhat
        Q0 % Preconditioner for G, size(Q)=[1 1 Nj x 1]. G*eta = etahat, (QG)*eta = Q*etahat, so etabar==Q*etahat.


        K2unique     % unique squared-wavenumbers
        iK2unique    % map from 2-dim K2, to 1-dim K2unique
        K2uniqueK2Map % cell array Nk in length. Each cell contains indices back to K2

        % IGW transformation matrices
        PFpmInv, QGpmInv % size(PFinv,PGinv)=[Nz x Nj x Nk]
        PFpm, QGpm % size(PF,PG)=[Nj x Nz x Nk]
        QGwg  % size(PF,PG)=[Nj x Nj x Nk]
        Ppm % Preconditioner for F, size(P)=[Nj x Nk]. F*u = uhat, (PF)*u = P*uhat, so ubar==P*uhat
        Qpm % Preconditioner for G, size(Q)=[Nj x Nk]. G*eta = etahat, (QG)*eta = Q*etahat, so etabar==Q*etahat.

        zInterp
        PFinvInterp, QGinvInterp

        Wzkl
        Wklz
        Aklz
    end

    properties (Dependent, SetAccess=private)
        nK2unique    % number of unique squared-wavenumbers
    end

    properties (GetAccess=public)
        h_pm % [Nj Nkl]
        h_0 % [Nj 1]
        iOmega
        isHydrostatic = 0
    end

    properties (Dependent)
        FinvMatrix
        GinvMatrix
        FMatrix
        GMatrix
    end
        
    methods
         
        function self = WVTransformBoussinesq(Lxyz, Nxyz, options)
            arguments
                Lxyz (1,3) double {mustBePositive}
                Nxyz (1,3) double {mustBePositive}
                options.rho function_handle = @isempty
                options.N2 function_handle = @isempty
                options.dLnN2func function_handle = @isempty
                options.latitude (1,1) double = 33
                options.rho0 (1,1) double {mustBePositive} = 1025
                options.shouldAntialias double = 1
                options.jAliasingFraction double {mustBePositive(options.jAliasingFraction),mustBeLessThanOrEqual(options.jAliasingFraction,1)} = 2/3

                % ALL of these must be set for direct initialization to
                % avoid actually computing the modes.
                options.dLnN2 (:,1) double
                options.PF0inv
                options.QG0inv
                options.PF0
                options.QG0
                options.h_0
                options.P0
                options.Q0
                options.z 
                options.PFpmInv
                options.QGpmInv
                options.PFpm
                options.QGpm
                options.h_pm
                options.Ppm 
                options.Qpm 
                options.QGwg 
            end

            % First we need to initialize the WVStratifiedFlow.
            if isfield(options,'z')
                z=options.z;
            else
                z = WVStratifiedFlow.quadraturePointsForStratifiedFlow(Lxyz(3),Nxyz(3),rho=options.rho,N2=options.N2,latitude=options.latitude);
            end
            self@WVStratifiedFlow(Lxyz(3),z,rho=options.rho,N2=options.N2,dLnN2=options.dLnN2func,latitude=options.latitude)

            % if all of these things are set initially (presumably read
            % from file), then we can initialize without computing modes.
            canInitializeDirectly = all(isfield(options,{'N2','latitude','rho0','dLnN2','PF0inv','QG0inv','PF0','QG0','h_0','P0','Q0','z','PFpmInv','QGpmInv','PFpm','QGpm','h_pm','Ppm','Qpm','QGwg'}));

            if canInitializeDirectly
                fprintf('Initialize the WVTransformBoussinesq directly from matrices.\n');
                Nj = size(options.PF0,1);
            else
                nModes = Nxyz(3)-1;
                if options.shouldAntialias == 1
                    Nj = floor(options.jAliasingFraction*nModes);
                else
                    Nj = nModes;
                end
            end

            self@WVTransform(Lxyz, Nxyz(1:2), z, latitude=options.latitude,rho0=options.rho0,Nj=Nj,shouldAntialias=options.shouldAntialias);

            if canInitializeDirectly
                self.PF0inv = options.PF0inv;
                self.QG0inv = options.QG0inv;
                self.PF0 = options.PF0;
                self.QG0 = options.QG0;
                self.h_0 = options.h_0;
                self.P0 = options.P0;
                self.Q0 = options.Q0;

                self.PFpmInv = options.PFpmInv;
                self.QGpmInv = options.QGpmInv;
                self.PFpm = options.PFpm;
                self.QGpm = options.QGpm;
                self.h_pm = options.h_pm;
                self.Ppm = options.Ppm;
                self.Qpm = options.Qpm;
                self.QGwg = options.QGwg;

                % K2unique are the unique wavenumbers (sorted)
                % iK2unique is the same size as K2, but are the indices for
                % the K2unique matrix to recreate/map back to K2unique.
                Kh = self.Kh;
                K2 = reshape((Kh(1,:)).^2,[],1);
                [self.K2unique,~,self.iK2unique] = unique(K2);
                self.iK2unique = reshape(self.iK2unique,size(K2));
                self.K2uniqueK2Map = cell(length(self.K2unique),1);
                for iK=1:length(self.K2unique)
                    self.K2uniqueK2Map{iK} = find(self.iK2unique==iK);
                end
            else
                % K2unique are the unique wavenumbers (sorted)
                % iK2unique is the same size as K2, but are the indices for
                % the K2unique matrix to recreate/map back to K2unique.
                Kh = self.Kh;
                K2 = reshape((Kh(1,:)).^2,[],1);
                [self.K2unique,~,self.iK2unique] = unique(K2);
                self.iK2unique = reshape(self.iK2unique,size(K2));
                self.K2uniqueK2Map = cell(length(self.K2unique),1);
                for iK=1:length(self.K2unique)
                    self.K2uniqueK2Map{iK} = find(self.iK2unique==iK);
                end

                self.buildVerticalModeProjectionOperators();
            end

            self.initializeStratifiedFlow();
            self.initializeGeostrophicComponent();
            self.initializeMeanDensityAnomalyComponent();
            self.initializeInternalGravityWaveComponent();
            self.initializeInertialOscillationComponent();

            % self.offgridModes = WVOffGridTransform(im,self.latitude, self.N2Function,1);

            self.addPropertyAnnotations(WVPropertyAnnotation('PF0inv',{'z','j'},'','Preconditioned F-mode inverse transformation'));
            self.addPropertyAnnotations(WVPropertyAnnotation('QG0inv',{'z','j'},'','Preconditioned G-mode inverse transformation'));
            self.addPropertyAnnotations(WVPropertyAnnotation('PF0',{'j','z'},'','Preconditioned F-mode forward transformation'));
            self.addPropertyAnnotations(WVPropertyAnnotation('QG0',{'j','z'},'','Preconditioned G-mode forward transformation'));
            self.addPropertyAnnotations(WVPropertyAnnotation('P0',{'j'},'','Preconditioner for F, size(P)=[1 Nj]. F*u = uhat, (PF)*u = P*uhat, so ubar==P*uhat'));
            self.addPropertyAnnotations(WVPropertyAnnotation('Q0',{'j'},'','Preconditioner for G, size(Q)=[1 Nj]. G*eta = etahat, (QG)*eta = Q*etahat, so etabar==Q*etahat. '));
            self.addPropertyAnnotations(WVPropertyAnnotation('h_0',{'j'},'m', 'equivalent depth of each geostrophic mode', detailedDescription='- topic: Domain Attributes — Stratification'));

            self.addDimensionAnnotations(WVDimensionAnnotation('K2unique', 'rad/m', 'unique horizontal wavenumbers (sorted)'));
            self.addPropertyAnnotations(WVPropertyAnnotation('iK2unique',{'kl'},'index', 'index for the K2 unique matrix'));

            self.addPropertyAnnotations(WVPropertyAnnotation('PFpmInv',{'z','j','K2unique'},'','Preconditioned F-mode inverse transformation'));
            self.addPropertyAnnotations(WVPropertyAnnotation('QGpmInv',{'z','j','K2unique'},'','Preconditioned G-mode inverse transformation'));
            self.addPropertyAnnotations(WVPropertyAnnotation('PFpm',{'j','z','K2unique'},'','Preconditioned F-mode forward transformation'));
            self.addPropertyAnnotations(WVPropertyAnnotation('QGpm',{'j','z','K2unique'},'','Preconditioned G-mode forward transformation'));
            self.addPropertyAnnotations(WVPropertyAnnotation('Ppm',{'j','K2unique'},'','Preconditioner for F, size(P)=[1 Nj]. F*u = uhat, (PF)*u = P*uhat, so ubar==P*uhat'));
            self.addPropertyAnnotations(WVPropertyAnnotation('Qpm',{'j','K2unique'},'','Preconditioner for G, size(Q)=[1 Nj]. G*eta = etahat, (QG)*eta = Q*etahat, so etabar==Q*etahat. '));
            self.addPropertyAnnotations(WVPropertyAnnotation('QGwg',{'j','j','K2unique'},'','Transformation from geostrophic to wave-modes'));
            self.addPropertyAnnotations(WVPropertyAnnotation('h_pm',{'j','kl'},'m', 'equivalent depth of each wave mode', detailedDescription='- topic: Domain Attributes — Stratification'));

            self.Aklz = zeros(self.Nx*self.Ny,self.Nz);
            self.Wklz = complex(zeros(self.Nx*self.Ny,self.Nz));
            self.Wzkl = complex(zeros(self.Nz,self.Nkl));

            self.nonlinearFluxOperation = WVNonlinearFlux(self);
        end

        function value = get.nK2unique(self)
            value=length(self.K2unique);
        end

        function wvtX2 = waveVortexTransformWithResolution(self,m)
            if ~isempty(self.dLnN2Function)
                wvtX2 = WVTransformBoussinesq([self.Lx self.Ly self.Lz],m, self.rhoFunction,latitude=self.latitude,rho0=self.rho0, N2func=self.N2Function, dLnN2func=self.dLnN2Function);
            else
                wvtX2 = WVTransformBoussinesq([self.Lx self.Ly self.Lz],m,latitude=self.latitude,rho0=self.rho0, N2=self.N2Function);
            end

            wvtX2.t0 = self.t0;
            [wvtX2.Ap,wvtX2.Am,wvtX2.A0] = self.spectralVariableWithResolution(wvtX2,self.Ap,self.Am,self.A0);
            wvtX2.nonlinearFluxOperation = self.nonlinearFluxOperation.nonlinearFluxWithResolutionOfTransform(wvtX2);
        end

        function self = buildVerticalModeProjectionOperators(self)
            [self.P0,self.Q0,self.PF0inv,self.PF0,self.QG0inv,self.QG0,self.h_0] = self.verticalProjectionOperatorsForGeostrophicModes(self.Nj);

            self.PFpmInv = zeros(self.Nz,self.Nj,self.nK2unique);
            self.QGpmInv = zeros(self.Nz,self.Nj,self.nK2unique);
            self.PFpm =    zeros(self.Nj,self.Nz,self.nK2unique);
            self.QGpm =    zeros(self.Nj,self.Nz,self.nK2unique);
            h = zeros(self.Nj,self.nK2unique);
            self.Ppm =     zeros(self.Nj,self.nK2unique);
            self.Qpm =     zeros(self.Nj,self.nK2unique);
            self.QGwg =    zeros(self.Nj,self.Nj,self.nK2unique);
            for iK=1:self.nK2unique
                [self.Ppm(:,iK),self.Qpm(:,iK),self.PFpmInv(:,:,iK),self.PFpm(:,:,iK),self.QGpmInv(:,:,iK),self.QGpm(:,:,iK),h(:,iK)] = self.verticalProjectionOperatorsForIGWModes(sqrt(self.K2unique(iK)),self.Nj);
                self.QGwg(:,:,iK ) = self.QGpm(:,:,iK )*self.QG0inv;
            end

            self.h_pm = zeros(self.spectralMatrixSize);
            for iK=1:size(self.h_pm,2)
                self.h_pm(:,iK) = h(:,self.iK2unique(iK));
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Transformations FROM the spatial domain
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

        function u_bar = transformFromSpatialDomainWithFio(self,u)
            u_bar = (self.PFpm(:,:,1 )*u)./self.Ppm(:,1);
        end

        function u_bar = transformFromSpatialDomainWithFg(self, u)
            u_bar = (self.PF0*u)./self.P0;
        end

        function w_bar = transformFromSpatialDomainWithGg(self, w)
            w_bar = (self.QG0*w)./self.Q0;
        end

        function w_bar = transformWithG_wg(self, w_bar )
            w_bar = self.Q0 .* w_bar;
            for iK=1:length(self.K2unique)
                indices = self.K2uniqueK2Map{iK};
                w_bar(:,indices) = (self.QGwg(:,:,iK) * w_bar(:,indices))./self.Qpm(:,iK);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Transformations to and from the spatial domain
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
                if ~isscalar(options.Apm) && ~isscalar(options.A0)
                    self.wvBuffer = self.PF0inv*(self.P0 .* options.A0);
                    for iK=1:length(self.K2unique)
                        indices = self.K2uniqueK2Map{iK};
                        self.wvBuffer(:,indices) = self.wvBuffer(:,indices) + self.PFpmInv(:,:,iK )*(self.Ppm(:,iK) .* options.Apm(:,indices));
                    end
                elseif ~isscalar(options.Apm)
                    for iK=1:length(self.K2unique)
                        indices = self.K2uniqueK2Map{iK};
                        self.wvBuffer(:,indices) = self.PFpmInv(:,:,iK )*(self.Ppm(:,iK) .* options.Apm(:,indices));
                    end
                else
                    self.wvBuffer = self.PF0inv*(self.P0 .* options.A0);
                end
                u = self.transformToSpatialDomainWithFourier(self.wvBuffer);  
            end
        end

        function w = transformToSpatialDomainWithG(self, options)
            arguments
                self WVTransform {mustBeNonempty}
                options.Apm double = 0
                options.A0 double = 0
            end

            if isscalar(options.Apm) && isscalar(options.A0)
                w = zeros(self.spatialMatrixSize);
            else
                if ~isscalar(options.Apm) && ~isscalar(options.A0)
                    self.wvBuffer = self.QG0inv*(self.Q0 .* options.A0);
                    for iK=1:length(self.K2unique)
                        indices = self.K2uniqueK2Map{iK};
                        self.wvBuffer(:,indices) = self.wvBuffer(:,indices) + self.QGpmInv(:,:,iK )*(self.Qpm(:,iK) .* options.Apm(:,indices));
                    end
                elseif ~isscalar(options.Apm)
                    for iK=1:length(self.K2unique)
                        indices = self.K2uniqueK2Map{iK};
                        self.wvBuffer(:,indices) = self.QGpmInv(:,:,iK )*(self.Qpm(:,iK) .* options.Apm(:,indices));
                    end
                else
                    self.wvBuffer = self.QG0inv*(self.Q0 .* options.A0);
                end
                w = self.transformToSpatialDomainWithFourier(self.wvBuffer);
            end
        end


        function u = transformToSpatialDomainWithFg(self, u_bar)
            % arguments
            %     self WVTransform {mustBeNonempty}
            %     u_bar
            % end
            u = self.PF0inv*(self.P0 .* u_bar);
        end

        function w = transformToSpatialDomainWithGg(self, w_bar)
            % arguments
            %     self WVTransform {mustBeNonempty}
            %     w_bar
            % end
            % simply changing QG0inv to PF0inv dramatically increases the
            % speed of the downstream function transformFromWVGridToDFTGrid.
            % Why?
            w = self.QG0inv*(self.Q0 .* w_bar);
        end
        
        function u = transformToSpatialDomainWithFw(self, u_bar)
            u = zeros(self.Nz,self.Nkl);
            for iK=1:length(self.K2unique)
                indices = self.K2uniqueK2Map{iK};
                u_bar(:,indices) = self.Ppm(:,iK) .* u_bar(:,indices);
                u(:,indices) = self.PFpmInv(:,:,iK )*u_bar(:,indices);
            end
        end
                
        function w = transformToSpatialDomainWithGw(self, w_bar)
            w = zeros(self.Nz,self.Nkl);
            for iK=1:length(self.K2unique)
                indices = self.K2uniqueK2Map{iK};
                w_bar(:,indices) = self.Qpm(:,iK) .* w_bar(:,indices);
                w(:,indices) = self.QGpmInv(:,:,iK )*w_bar(:,indices);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Higher resolution vertical grid for interpolation
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = buildInterpolationProjectionOperators(self,dof)
            zInterp_ = cat(1,-self.Lz,-self.Lz + cumsum(reshape(shiftdim(repmat(diff(self.z)/dof,[1 dof]),1),[],1)));
            zInterp_(end) = self.z(end);
            self.buildInterpolationProjectionOperatorsForGrid(zInterp_);
        end

        function self = buildInterpolationProjectionOperatorsForGrid(self,zInterp)
            self.zInterp = zInterp;
            im = InternalModesWKBSpectral(N2=self.N2Function,zIn=[-self.Lz 0],zOut=self.zInterp,latitude=self.latitude,nModes=self.verticalModes.nModes);
            im.normalization = Normalization.geostrophicFreeSurface;
            im.upperBoundary = UpperBoundary.rigidLid;
            [Finv,Ginv] = im.ModesAtFrequency(0);
            N = length(self.zInterp);

            % dump the Nyquist mode
            Finv = Finv(:,1:end-1);
            Ginv = Ginv(:,1:end-1);

            % add the barotropic mode
            Finv = cat(2,ones(N,1),Finv);
            Ginv = cat(2,zeros(N,1),Ginv);

            self.PFinvInterp = Finv./shiftdim(self.P0,1);
            self.QGinvInterp = Ginv./shiftdim(self.Q0,1);

            self.addDimensionAnnotations(WVDimensionAnnotation('z-interp', 'm', 'z-coordinate dimension for interpolation'));

            outputVar = WVVariableAnnotation('uInterp',{'x','y','z-interp'},'m/s', 'x-component of the fluid velocity');
            f = @(wvt) wvt.transformToSpatialDomainWithFInterp(wvt.UAp.*wvt.Apt + wvt.UAm.*wvt.Amt + wvt.UA0.*wvt.A0t);
            self.addOperation(WVOperation('uInterp',outputVar,f));

            outputVar = WVVariableAnnotation('vInterp',{'x','y','z-interp'},'m/s', 'y-component of the fluid velocity');
            f = @(wvt) wvt.transformToSpatialDomainWithFInterp(wvt.VAp.*wvt.Apt + wvt.VAm.*wvt.Amt + wvt.VA0.*wvt.A0t);
            self.addOperation(WVOperation('vInterp',outputVar,f));

            outputVar = WVVariableAnnotation('wInterp',{'x','y','z-interp'},'m/s', 'z-component of the fluid velocity');
            f = @(wvt) wvt.transformToSpatialDomainWithGInterp(wvt.WAp.*wvt.Apt + wvt.WAm.*wvt.Amt);
            self.addOperation(WVOperation('wInterp',outputVar,f));

            outputVar = WVVariableAnnotation('pInterp',{'x','y','z-interp'},'kg/m/s2', 'pressure anomaly');
            f = @(wvt) wvt.rho0*wvt.g*wvt.transformToSpatialDomainWithFInterp(wvt.NAp.*wvt.Apt + wvt.NAm.*wvt.Amt + wvt.NA0.*wvt.A0t);
            self.addOperation(WVOperation('pInterp',outputVar,f));

            outputVar = WVVariableAnnotation('etaInterp',{'x','y','z-interp'},'m', 'isopycnal deviation');
            f = @(wvt) wvt.transformToSpatialDomainWithGInterp(wvt.NAp.*wvt.Apt + wvt.NAm.*wvt.Amt + wvt.NA0.*wvt.A0t);
            self.addOperation(WVOperation('etaInterp',outputVar,f));
        end

        function u = transformToSpatialDomainWithFInterp(self, u_bar)
            u_bar = (self.P0 .* u_bar)*(self.Nx*self.Ny);
            % hydrostatic modes commute with the DFT
            u_bar = ifft(ifft(u_bar,self.Nx,1),self.Ny,2,'symmetric');
            u_bar = permute(u_bar,[3 1 2]); % keep adjacent in memory
            u_bar = reshape(u_bar,self.Nj,[]);
            u = self.PFinvInterp*u_bar;
            u = reshape(u,length(self.zInterp),self.Nx,self.Ny);
            u = permute(u,[2 3 1]);
        end

        function w = transformToSpatialDomainWithGInterp(self, w_bar )
            w_bar = (self.Q0 .* w_bar)*(self.Nx*self.Ny);
            % hydrostatic modes commute with the DFT
            w_bar = ifft(ifft(w_bar,self.Nx,1),self.Ny,2,'symmetric');

            w_bar = permute(w_bar,[3 1 2]); % keep adjacent in memory
            w_bar = reshape(w_bar,self.Nj,[]);
            w = self.QGinvInterp*w_bar;
            w = reshape(w,length(self.zInterp),self.Nx,self.Ny);
            w = permute(w,[2 3 1]);
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

            Finv = shiftdim(wvt.P0,1) .* wvt.PF0inv;

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

            F = wvt.PF0 ./ shiftdim(wvt.P0,2);

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

            Ginv = shiftdim(wvt.Q0,1) .* wvt.QG0inv;

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

            G = wvt.QG0 ./ shiftdim(wvt.Q0,2);

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
            ratio = self.P0(j+1);
        end

        function ratio = maxFg(self,kMode,lMode,j)
            arguments
                self WVTransform {mustBeNonempty}
                kMode (:,1) double
                lMode (:,1) double
                j (:,1) double
            end
            ratio = self.P0(j+1);
        end

        [ncfile,matFilePath] = writeToFile(wvt,path,variables,options)

        function flag = isequal(self,other)
            arguments
                self WVTransform
                other WVTransform
            end
            flag = isequal@WVTransform(self,other);
            flag = flag & isequal(self.dLnN2, other.dLnN2);
            flag = flag & isequal(self.PF0inv, other.PF0inv);
            flag = flag & isequal(self.QG0inv, other.QG0inv);
            flag = flag & isequal(self.PF0,other.PF0);
            flag = flag & isequal(self.QG0,other.QG0);
            flag = flag & isequal(self.P0, other.P0);
            flag = flag & isequal(self.Q0, other.Q0);
            flag = flag & isequal(self.h_0, other.h_0);
            flag = flag & isequal(self.PFpmInv, other.PFpmInv);
            flag = flag & isequal(self.QGpmInv, other.QGpmInv);
            flag = flag & isequal(self.PFpm,other.PFpm);
            flag = flag & isequal(self.QGpm,other.QGpm);
            flag = flag & isequal(self.Ppm, other.Ppm);
            flag = flag & isequal(self.Qpm, other.Qpm);
            flag = flag & isequal(self.h_pm, other.h_pm);
            flag = flag & isequal(self.QGwg, other.QGwg);
        end
    end

    methods (Static)
        wvt = waveVortexTransformFromFile(path,options)
    end
   
        
        
end 



