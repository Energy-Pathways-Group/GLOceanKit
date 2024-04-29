classdef WVTransformBoussinesq < WVTransform
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

    properties (GetAccess=public, SetAccess=protected)
        rhobar, N2, dLnN2 % on the z-grid, size(N2) = [length(z) 1];
        rhoFunction, N2Function, dLnN2Function % function handles

        internalModes

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
        
    methods
         
        function self = WVTransformBoussinesq(Lxyz, Nxyz, options)
            arguments
                Lxyz (1,3) double {mustBePositive}
                Nxyz (1,3) double {mustBePositive}
                options.rho function_handle
                options.N2 function_handle
                options.dLnN2func function_handle = @disp
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


            % if all of these things are set initially (presumably read
            % from file), then we can initialize without computing modes.
            canInitializeDirectly = all(isfield(options,{'N2','latitude','rho0','dLnN2','PF0inv','QG0inv','PF0','QG0','h_0','P0','Q0','z','PFpmInv','QGpmInv','PFpm','QGpm','h_pm','Ppm','Qpm','QGwg'}));

            if canInitializeDirectly
                fprintf('Initialize the WVTransformHydrostatic directly from matrices.\n');
                % We already know the quadrature points
                z = options.z;
                N2 = options.N2(z);
                im = InternalModesWKBSpectral(N2=options.N2,zIn=[-Lxyz(3) 0],zOut=z,latitude=options.latitude,nModes=Nxyz(3)-1);
                Nj = size(options.PF0,1);
            else
                nModes = Nxyz(3)-1;
                % Before initializing the superclass, we need to find the
                % Gauss-quadrature points for this stratification profile.
                Nz = Nxyz(3);
                z = linspace(-Lxyz(3),0,Nz*10)';
                if isfield(options,'N2')
                    im = InternalModesWKBSpectral(N2=options.N2,zIn=[-Lxyz(3) 0],zOut=z,latitude=options.latitude, nEVP=max(256,floor(2.1*Nz)));
                elseif isfield(options,'rho')
                    im = InternalModesWKBSpectral(rho=options.rho,zIn=[-Lxyz(3) 0],zOut=z,latitude=options.latitude);
                else
                    error('You must specify either rho or N2.');
                end
                im.normalization = Normalization.kConstant;
                im.upperBoundary = UpperBoundary.rigidLid;
                z = im.GaussQuadraturePointsForModesAtFrequency(nModes+1,im.f0);

                % If the user requests nModes---we will have that many fully
                % resolved modes. It works as follows for rigid lid:
                % - There is one barotropic mode that appears in F
                % - There are nModes-1 *internal modes* for G and F.
                % - We compute the nModes+1 internal mode for F, to make it
                % complete.
                % This is nModes+1 grid points necessary to make this happen.
                % This should make sense because there are nModes-1 internal
                % modes, but the boundaries.
                if isfield(options,'N2')
                    im = InternalModesWKBSpectral(N2=options.N2,zIn=[-Lxyz(3) 0],zOut=z,latitude=options.latitude,nModes=nModes,nEVP=128);
                    N2 = options.N2(z);
                    N2func = options.N2;
                    rhoFunc = im.rho_function;
                elseif isfield(options,'rho')
                    im = InternalModesWKBSpectral(rho=options.rho,zIn=[-Lxyz(3) 0],zOut=z,latitude=options.latitude,nModes=nModes,nEVP=128);
                    N2 = im.N2;
                    N2func = im.N2_function;
                    rhoFunc = options.rho;
                end

                if options.shouldAntialias == 1
                    Nj = floor(options.jAliasingFraction*nModes);
                else
                    Nj = nModes;
                end
            end
            im.normalization = Normalization.kConstant;
            im.upperBoundary = UpperBoundary.rigidLid;

            % This is enough information to initialize the superclass
            self@WVTransform(Lxyz, Nxyz(1:2), z, latitude=options.latitude,rho0=options.rho0,Nj=Nj,Nmax=sqrt(max(N2)));

            if canInitializeDirectly
                self.rhoFunction = im.rho_function;
                self.N2Function = options.N2;
                self.internalModes = im;

                self.rhobar = self.rhoFunction(self.z);
                self.N2 = self.N2Function(self.z);
                self.dLnN2 = options.dLnN2;

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

                self.buildTransformationMatrices();
            else
                if isequal(options.dLnN2func,@disp)
                    dLnN2 = im.rho_zz./im.rho_z;
                else
                    dLnN2 = options.dLnN2func(z);
                end

                self.rhoFunction = rhoFunc;
                self.N2Function = N2func;
                self.internalModes = im;

                self.rhobar = rhoFunc(self.z);
                self.N2 = N2;
                self.dLnN2 = dLnN2;
                
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

                self.BuildProjectionOperators();
            end

            self.offgridModes = WVOffGridTransform(im,self.latitude, self.N2Function,1);

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


            self.nonlinearFluxOperation = WVNonlinearFlux(self);
        end

        function value = get.nK2unique(self)
            value=length(self.K2unique);
        end

        function wvtX2 = waveVortexTransformWithResolution(self,m)
            if ~isempty(self.dLnN2Function)
                wvtX2 = WVTransformHydrostatic([self.Lx self.Ly self.Lz],m, self.rhoFunction,latitude=self.latitude,rho0=self.rho0, N2func=self.N2Function, dLnN2func=self.dLnN2Function);
            else
                wvtX2 = WVTransformHydrostatic([self.Lx self.Ly self.Lz],m,latitude=self.latitude,rho0=self.rho0, N2=self.N2Function);
            end

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

        function self = BuildProjectionOperators(self)
            self.PF0inv = zeros(self.Nz,self.Nj,1);
            self.QG0inv = zeros(self.Nz,self.Nj,1);
            self.PF0 =    zeros(self.Nj,self.Nz,1);
            self.QG0 =    zeros(self.Nj,self.Nz,1);
            self.h_0 =    zeros(self.Nj,1);
            self.P0 =     zeros(self.Nj,1);
            self.Q0 =     zeros(self.Nj,1);

            [self.P0,self.Q0,self.PF0inv,self.PF0,self.QG0inv,self.QG0,self.h_0] = self.BuildProjectionOperatorsForGeostrophicModes();

            self.PFpmInv = zeros(self.Nz,self.Nj,self.nK2unique);
            self.QGpmInv = zeros(self.Nz,self.Nj,self.nK2unique);
            self.PFpm =    zeros(self.Nj,self.Nz,self.nK2unique);
            self.QGpm =    zeros(self.Nj,self.Nz,self.nK2unique);
            h = zeros(self.Nj,self.nK2unique);
            self.Ppm =     zeros(self.Nj,self.nK2unique);
            self.Qpm =     zeros(self.Nj,self.nK2unique);
            self.QGwg =    zeros(self.Nj,self.Nj,self.nK2unique);
            for iK=1:self.nK2unique
                [self.Ppm(:,iK),self.Qpm(:,iK),self.PFpmInv(:,:,iK),self.PFpm(:,:,iK),self.QGpmInv(:,:,iK),self.QGpm(:,:,iK),h(:,iK)] = self.BuildProjectionOperatorsForIGWModes(sqrt(self.K2unique(iK)));
                self.QGwg(:,:,iK ) = self.QGpm(:,:,iK )*self.QG0inv;
            end

            self.h_pm = zeros(self.spectralMatrixSize);
            for iK=1:size(self.h_pm,2)
                self.h_pm(:,iK) = h(:,self.iK2unique(iK));
            end

            self.buildTransformationMatrices();
        end

        function [P,Q,PFinv,PF,QGinv,QG,h] = BuildProjectionOperatorsForGeostrophicModes(self)
            % Now go compute the appropriate number of modes at the
            % quadrature points.
            self.internalModes.normalization = Normalization.kConstant;
            self.internalModes.upperBoundary = UpperBoundary.rigidLid;
            [Finv,Ginv,h] = self.internalModes.ModesAtFrequency(0);
            [P,Q,PFinv,PF,QGinv,QG,h] = self.BuildProjectionOperatorsWithRigidLid(Finv,Ginv,h);
        end

        function [P,Q,PFinv,PF,QGinv,QG,h] = BuildProjectionOperatorsForIGWModes(self,k)
            % Now go compute the appropriate number of modes at the
            % quadrature points.
            self.internalModes.normalization = Normalization.kConstant;
            self.internalModes.upperBoundary = UpperBoundary.rigidLid;
            [Finv,Ginv,h] = self.internalModes.ModesAtWavenumber(k);
            [P,Q,PFinv,PF,QGinv,QG,h] = self.BuildProjectionOperatorsWithRigidLid(Finv,Ginv,h);
        end

        function [P,Q,PFinv,PF,QGinv,QG,h] = BuildProjectionOperatorsWithRigidLid(self,Finv,Ginv,h)
            nModes = size(Finv,2);

            % Make these matrices invertible by adding the barotropic mode
            % to F, and removing the boundaries of G.
            Finv = cat(2,ones(self.Nz,1),Finv);
            Ginv = Ginv(2:end-1,1:end-1);

            % Compute the precondition matrices (really, diagonals)
            P = max(abs(Finv),[],1); % ones(1,size(Finv,1)); %
            Q = max(abs(Ginv),[],1); % ones(1,size(Ginv,1)); %

            % Now create the actual transformation matrices
            PFinv = Finv./P;
            QGinv = Ginv./Q;
            PF = inv(PFinv);
            QG = inv(QGinv);

            maxCond = max([cond(PFinv), cond(QGinv), cond(PF), cond(QG)],[],2);
            if maxCond > 1000
                warning('Condition number is %f the vertical transformations.',maxCond);
            end
            % size(F)=[Nz x Nj+1], barotropic mode AND extra Nyquist mode
            % but, we will only multiply by vectors [Nj 1], so dump the
            % last column. Now size(Fp) = [Nz x Nj].
            PFinv = PFinv(:,1:end-1);

            % size(Finv)=[Nj+1, Nz], but we don't care about the last mode
            PF = PF(1:end-1,:);

            % size(G) = [Nz-2, Nj-1], need zeros for the boundaries
            % and add the 0 barotropic mode, so size(G) = [Nz, Nj],
            QGinv = cat(2,zeros(self.Nz,1),cat(1,zeros(1,nModes-1),QGinv,zeros(1,nModes-1)));

            % size(Ginv) = [Nj-1, Nz-2], need a zero for the barotropic
            % mode, but also need zeros for the boundary
            QG = cat(2,zeros(nModes,1), cat(1,zeros(1,self.Nz-2),QG),zeros(nModes,1));

            % want size(h)=[Nj 1]
            h = cat(1,1,reshape(h(1:end-1),[],1)); % remove the extra mode at the end

            P = reshape(P(1:end-1),[],1);
            Q = reshape(cat(2,1,Q),[],1);

            PFinv = PFinv(:,1:self.Nj);
            PF = PF(1:self.Nj,:);
            P = P(1:self.Nj,1);
            QGinv = QGinv(:,1:self.Nj);
            QG = QG(1:self.Nj,:);
            Q = Q(1:self.Nj,1);
            h = h(1:self.Nj,1);
        end

        function [P,Q,PFinv,PF,QGinv,QG,h] = BuildProjectionOperatorsWithFreeSurface(self,Finv,Ginv,h)
            % Make these matrices invertible by adding the barotropic mode
            % to F, and removing the lower boundary of G.
            Finv = cat(2,ones(self.Nz,1),Finv); % [Nz Nj+1]
            Ginv = Ginv(2:end,:);  % [Nz-1 Nj]

            % Compute the precondition matrices (really, diagonals)
            P = max(abs(Finv),[],1); % ones(1,size(Finv,1)); %
            Q = max(abs(Ginv),[],1); % ones(1,size(Ginv,1)); %

            % Now create the actual transformation matrices
            PFinv = Finv./P;
            QGinv = Ginv./Q;
            PF = inv(PFinv); % [Nj+1 Nz]
            QG = inv(QGinv); % [Nj Nz-1]

            maxCond = max([cond(PFinv), cond(QGinv), cond(PF), cond(QG)],[],2);
            if maxCond > 1000
                warning('Condition number is %f the vertical transformations.',maxCond);
            end
            % size(PFinv)=[Nz x Nj+1], barotropic mode AND extra Nyquist mode
            % but, we will only multiply by vectors [Nj 1], so dump the
            % last column. Now size(PFinv) = [Nz x Nj].
            PFinv = PFinv(:,1:end-1);

            % size(PF)=[Nj+1, Nz], but we don't care about the last mode
            PF = PF(1:end-1,:);

            % size(QGinv) = [Nz-1, Nj], need zeros for the lower boundaries
            % and add the 0 barotropic mode, so size(G) = [Nz, Nj],
            QGinv = QGinv(:,1:end-1); % dump Nyquist  
            QGinv = cat(2,zeros(self.Nz-1,1),QGinv); % add barotropic mode
            QGinv = cat(1,zeros(1,self.Nj),QGinv); % add zeros at along the bottom

            % Now have to do the same thing to the condition matrix
            Q = cat(2,0,Q(1:end-1));

            % size(QG) = [Nj, Nz-1], need a zero for the barotropic
            % mode, but also need zeros for the boundary
            QG = cat(1,zeros(1,self.Nz-1),QG(1:end-1,:)); % dump the Nyquist mode, add a barotropic mode (all zeros)
            QG = cat(2,zeros(self.Nj,1), QG); % add the bottom boundary

            % want size(h)=[1 1 Nj]
            h = shiftdim(h,-1);

            P = shiftdim(P(1:end-1),-1);
            Q = shiftdim(Q,-1);
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
                w_bar(:,indices) = self.QGwg(:,:,iK) * w_bar(:,indices);
                w_bar(:,indices) = w_bar(:,indices)./self.Qpm(:,iK);
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
                options.Apm double = []
                options.A0 double = []
            end
            u_bar = 0;
            if ~isempty(options.Apm)
                u_bar = u_bar + self.transformToSpatialDomainWithFw(options.Apm);
            end
            if ~isempty(options.A0)
                u_bar = u_bar + self.transformToSpatialDomainWithFg(options.A0);
            end
            u = self.transformToSpatialDomainWithFourier(self.horizontalGeometry.transformFromWVGridToDFTGrid(u_bar));
        end

        function u = transformToSpatialDomainWithG(self, options)
            arguments
                self WVTransform {mustBeNonempty}
                options.Apm double = []
                options.A0 double  = []
            end
            u_bar = 0;
            if ~isempty(options.Apm)
                u_bar = u_bar + self.transformToSpatialDomainWithGw(options.Apm);
            end
            if ~isempty(options.A0)
                u_bar = u_bar + self.transformToSpatialDomainWithGg(options.A0);
            end
            u = self.transformToSpatialDomainWithFourier(self.horizontalGeometry.transformFromWVGridToDFTGrid(u_bar));
        end

        function u = transformToSpatialDomainWithFg(self, u_bar)
            arguments
                self WVTransform {mustBeNonempty}
                u_bar
            end
            u = self.PF0inv*(self.P0 .* u_bar);
        end

        function w = transformToSpatialDomainWithGg(self, w_bar)
            arguments
                self WVTransform {mustBeNonempty}
                w_bar
            end
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
            im = InternalModesWKBSpectral(N2=self.N2Function,zIn=[-self.Lz 0],zOut=self.zInterp,latitude=self.latitude,nModes=self.internalModes.nModes);
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
        % Needed to add and remove internal waves from the model
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function ratio = uMaxGNormRatioForWave(self,k0, l0, j0)
            ratio = 1/self.P0(j0+1);
        end   

        function ratio = uMaxA0(self,k0, l0, j0)
            error('This is pulled straight from hydrostatic and must be updated.')
            if j0 == 0
                ratio = 1;
            else
                ratio = abs(1/self.F_g(k0+1,l0+1,j0+1));
            end
        end

        [ncfile,matFilePath] = writeToFile(wvt,path,variables,options)

    end

    methods (Static)
        wvt = waveVortexTransformFromFile(path,options)
    end
   
        
        
end 



