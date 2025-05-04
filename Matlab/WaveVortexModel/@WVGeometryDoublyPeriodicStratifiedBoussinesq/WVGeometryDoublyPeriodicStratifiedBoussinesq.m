classdef WVGeometryDoublyPeriodicStratifiedBoussinesq < WVGeometryDoublyPeriodicStratified
    properties (Access=public) %(GetAccess=private, SetAccess=private) %(Access=private)
        K2unique     % unique squared-wavenumbers
        iK2unique    % map from 2-dim K2, to 1-dim K2unique
        K2uniqueK2Map % cell array Nk in length. Each cell contains indices back to K2

        % IGW transformation matrices
        PFpmInv, QGpmInv % size(PFinv,PGinv)=[Nz x Nj x Nk]
        PFpm, QGpm % size(PF,PG)=[Nj x Nz x Nk]
        QGwg  % size(PF,PG)=[Nj x Nj x Nk]
        Ppm % Preconditioner for F, size(P)=[Nj x Nk]. F*u = uhat, (PF)*u = P*uhat, so ubar==P*uhat
        Qpm % Preconditioner for G, size(Q)=[Nj x Nk]. G*eta = etahat, (QG)*eta = Q*etahat, so etabar==Q*etahat.

        wvBuffer
    end
    
    properties (Dependent, SetAccess=private)
        nK2unique    % number of unique squared-wavenumbers
    end

    methods
        function self = WVGeometryDoublyPeriodicStratifiedBoussinesq(Lxyz, Nxyz, options, directInit)
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

                options.shouldAntialias (1,1) logical = true
                options.z (:,1) double {mustBeNonempty} % quadrature points!
                options.j (:,1) double {mustBeNonempty}
                options.Nj (1,1) double {mustBePositive}
                options.rhoFunction function_handle = @isempty
                options.N2Function function_handle = @isempty
                options.rho0 (1,1) double {mustBePositive} = 1025
                options.planetaryRadius (1,1) double = 6.371e6
                options.rotationRate (1,1) double = 7.2921E-5
                options.latitude (1,1) double = 33
                options.g (1,1) double = 9.81
                options.dLnN2 (:,1) double
                options.PF0inv
                options.QG0inv
                options.PF0
                options.QG0
                options.h_0 (:,1) double
                options.P0 (:,1) double
                options.Q0 (:,1) double
                options.z_int (:,1) double

                directInit.PFpmInv
                directInit.QGpmInv
                directInit.PFpm
                directInit.QGpm
                directInit.h_pm
                directInit.Ppm
                directInit.Qpm
                directInit.QGwg
                directInit.K2unique
                directInit.iK2unique
            end

            optionCell = namedargs2cell(options);
            self@WVGeometryDoublyPeriodicStratified(Lxyz, Nxyz, optionCell{:})

            allFields = cell2struct([struct2cell(options);struct2cell(directInit)],[fieldnames(options);fieldnames(directInit)]);
            canInitializeDirectly = all(isfield(allFields, WVGeometryDoublyPeriodicStratifiedBoussinesq.newRequiredPropertyNames));

            if canInitializeDirectly
                self.PFpmInv = directInit.PFpmInv;
                self.QGpmInv = directInit.QGpmInv;
                self.PFpm = directInit.PFpm;
                self.QGpm = directInit.QGpm;
                self.h_pm = directInit.h_pm;
                self.Ppm = directInit.Ppm;
                self.Qpm = directInit.Qpm;
                self.QGwg = directInit.QGwg;
                self.K2unique = directInit.K2unique;
                self.iK2unique = directInit.iK2unique;

                % K2unique are the unique wavenumbers (sorted)
                % iK2unique is the same size as K2, but are the indices for
                % the K2unique matrix to recreate/map back to K2unique.
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

            self.wvBuffer = zeros(self.Nz,self.Nkl);
        end

        function self = buildVerticalModeProjectionOperators(self)
            self.PFpmInv = zeros(self.Nz,self.Nj,self.nK2unique);
            self.QGpmInv = zeros(self.Nz,self.Nj,self.nK2unique);
            self.PFpm =    zeros(self.Nj,self.Nz,self.nK2unique);
            self.QGpm =    zeros(self.Nj,self.Nz,self.nK2unique);
            h = zeros(self.Nj,self.nK2unique);
            self.Ppm =     zeros(self.Nj,self.nK2unique);
            self.Qpm =     zeros(self.Nj,self.nK2unique);
            self.QGwg =    zeros(self.Nj,self.Nj,self.nK2unique);

            fprintf("building %d boussinesq transformation matrices...",self.nK2unique)
            for iK=1:self.nK2unique
                if mod(iK,50) == 1
                    fprintf("%d...",iK);
                end
                [self.Ppm(:,iK),self.Qpm(:,iK),self.PFpmInv(:,:,iK),self.PFpm(:,:,iK),self.QGpmInv(:,:,iK),self.QGpm(:,:,iK),h(:,iK)] = self.verticalProjectionOperatorsForIGWModes(sqrt(self.K2unique(iK)),self.Nj);
                self.QGwg(:,:,iK ) = self.QGpm(:,:,iK )*self.QG0inv;
            end
            fprintf("finished.\n");

            self.h_pm = zeros(self.spectralMatrixSize);
            for iK=1:size(self.h_pm,2)
                self.h_pm(:,iK) = h(:,self.iK2unique(iK));
            end
        end

        function value = get.nK2unique(self)
            value=length(self.K2unique);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Transformation matrices
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

            [kMode,lMode] = self.primaryKLModeNumberFromKLModeNumber(kMode,lMode);
            iK = self.indexFromKLModeNumber(kMode,lMode);
            for iUnique=1:length(self.K2unique)
                if ismember(iK, self.K2uniqueK2Map{iUnique})
                    break
                end
            end
            Finv = shiftdim(self.Ppm(:,iUnique),1) .* self.PFpmInv(:,:,iUnique ); % 
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

            [kMode,lMode] = self.primaryKLModeNumberFromKLModeNumber(kMode,lMode);
            iK = self.indexFromKLModeNumber(kMode,lMode);
            for iUnique=1:length(self.K2unique)
                if ismember(iK, self.K2uniqueK2Map{iUnique})
                    break
                end
            end
            F = self.PFpm(:,:,iUnique )./self.Ppm(:,iUnique);
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

            [kMode,lMode] = self.primaryKLModeNumberFromKLModeNumber(kMode,lMode);
            iK = self.indexFromKLModeNumber(kMode,lMode);
            for iUnique=1:length(self.K2unique)
                if ismember(iK, self.K2uniqueK2Map{iUnique})
                    break
                end
            end
            Ginv = shiftdim(self.Qpm(:,iUnique),1) .* self.QGpmInv(:,:,iUnique ); % 
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

            [kMode,lMode] = self.primaryKLModeNumberFromKLModeNumber(kMode,lMode);
            iK = self.indexFromKLModeNumber(kMode,lMode);
            for iUnique=1:length(self.K2unique)
                if ismember(iK, self.K2uniqueK2Map{iUnique})
                    break
                end
            end
            G = self.QGpm(:,:,iUnique )./self.Qpm(:,iUnique);
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
            [kMode,lMode] = self.primaryKLModeNumberFromKLModeNumber(kMode,lMode);
            iK = self.indexFromKLModeNumber(kMode,lMode);
            for iUnique=1:length(self.K2unique)
                if ismember(iK, self.K2uniqueK2Map{iUnique})
                    break
                end
            end
            ratio = self.Ppm(j+1,iUnique);
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

        function w_bar = transformFromSpatialDomainWithG_w(self, w_hat )
            w_bar = zeros(self.spectralMatrixSize);
            for iK=1:length(self.K2unique)
                indices = self.K2uniqueK2Map{iK};
                w_bar(:,indices) = (self.QGpm(:,:,iK) * w_hat(:,indices))./self.Qpm(:,iK);
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
    end

    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % CAAnnotatedClass required methods, which enables writeToFile
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = WVGeometryDoublyPeriodicStratifiedBoussinesq.propertyAnnotationsForGeometry();
        end

        function vars = classRequiredPropertyNames()
            vars = WVGeometryDoublyPeriodicStratifiedBoussinesq.namesOfRequiredPropertiesForGeometry();
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Stratification specific property annotations and initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function requiredPropertyNames = namesOfRequiredPropertiesForGeometry()
            requiredPropertyNames = WVGeometryDoublyPeriodicStratified.namesOfRequiredPropertiesForGeometry();
            requiredPropertyNames = union(requiredPropertyNames,WVGeometryDoublyPeriodicStratifiedBoussinesq.newRequiredPropertyNames());
        end

        function newRequiredPropertyNames = newRequiredPropertyNames()
            newRequiredPropertyNames = {'K2unique','iK2unique','PFpmInv','QGpmInv','PFpm','QGpm','Ppm','Qpm','QGwg','h_pm'};
        end

        function propertyAnnotations = propertyAnnotationsForGeometry()
            propertyAnnotations = WVGeometryDoublyPeriodicStratified.propertyAnnotationsForGeometry();

            propertyAnnotations(end+1) = CADimensionProperty('K2unique', 'rad/m', 'unique horizontal wavenumbers (sorted)');
            propertyAnnotations(end+1) = CANumericProperty('iK2unique',{'kl'},'index', 'index for the K2 unique matrix');

            propertyAnnotations(end+1) = CANumericProperty('PFpmInv',{'z','j','K2unique'},'','Preconditioned F-mode inverse transformation');
            propertyAnnotations(end+1) = CANumericProperty('QGpmInv',{'z','j','K2unique'},'','Preconditioned G-mode inverse transformation');
            propertyAnnotations(end+1) = CANumericProperty('PFpm',{'j','z','K2unique'},'','Preconditioned F-mode forward transformation');
            propertyAnnotations(end+1) = CANumericProperty('QGpm',{'j','z','K2unique'},'','Preconditioned G-mode forward transformation');
            propertyAnnotations(end+1) = CANumericProperty('Ppm',{'j','K2unique'},'','Preconditioner for F, size(P)=[1 Nj]. F*u = uhat, (PF)*u = P*uhat, so ubar==P*uhat');
            propertyAnnotations(end+1) = CANumericProperty('Qpm',{'j','K2unique'},'','Preconditioner for G, size(Q)=[1 Nj]. G*eta = etahat, (QG)*eta = Q*etahat, so etabar==Q*etahat. ');
            propertyAnnotations(end+1) = CANumericProperty('QGwg',{'j','j','K2unique'},'','Transformation from geostrophic to wave-modes');
            propertyAnnotations(end+1) = CANumericProperty('h_pm',{'j','kl'},'m', 'equivalent depth of each wave mode', detailedDescription='- topic: Domain Attributes — Stratification');
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
            [Lxyz, Nxyz, geomOptions] = WVGeometryDoublyPeriodicStratified.requiredPropertiesForGeometryFromGroup(group);
            vars = CAAnnotatedClass.propertyValuesFromGroup(group,WVGeometryDoublyPeriodicStratifiedBoussinesq.newRequiredPropertyNames);
            newOptions = namedargs2cell(vars);
            options = cat(2,geomOptions,newOptions);
        end

        function geometry = geometryFromFile(path)
            arguments (Input)
                path char {mustBeNonempty}
            end
            arguments (Output)
                geometry WVGeometryDoublyPeriodicStratifiedBoussinesq {mustBeNonempty}
            end
            ncfile = NetCDFFile(path);
            geometry = WVGeometryDoublyPeriodicStratifiedBoussinesq.geometryFromGroup(ncfile);
        end

        function geometry = geometryFromGroup(group)
            arguments (Input)
                group NetCDFGroup {mustBeNonempty}
            end
            arguments (Output)
                geometry WVGeometryDoublyPeriodicStratifiedBoussinesq {mustBeNonempty}
            end
            CAAnnotatedClass.throwErrorIfMissingProperties(group,WVGeometryDoublyPeriodicStratifiedBoussinesq.namesOfRequiredPropertiesForGeometry);
            [Lxyz, Nxyz, options] = WVGeometryDoublyPeriodicStratifiedBoussinesq.requiredPropertiesForGeometryFromGroup(group);
            geometry = WVGeometryDoublyPeriodicStratifiedBoussinesq(Lxyz,Nxyz,options{:});
        end

    end
end