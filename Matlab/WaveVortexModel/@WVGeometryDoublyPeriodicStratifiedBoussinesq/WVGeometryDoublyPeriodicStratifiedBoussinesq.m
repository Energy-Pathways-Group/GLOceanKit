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
            end

            optionCell = namedargs2cell(options);
            self@WVGeometryDoublyPeriodicStratified(Lxyz, Nxyz, optionCell{:})

            allFields = cell2struct([struct2cell(options);struct2cell(directInit)],[fieldnames(options);fieldnames(directInit)]);
            canInitializeDirectly = all(isfield(allFields, WVStratification.namesOfRequiredPropertiesForStratification));

            if canInitializeDirectly
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

            for iK=1:self.nK2unique
                [self.Ppm(:,iK),self.Qpm(:,iK),self.PFpmInv(:,:,iK),self.PFpm(:,:,iK),self.QGpmInv(:,:,iK),self.QGpm(:,:,iK),h(:,iK)] = self.verticalProjectionOperatorsForIGWModes(sqrt(self.K2unique(iK)),self.Nj);
                self.QGwg(:,:,iK ) = self.QGpm(:,:,iK )*self.QG0inv;
            end

            self.h_pm = zeros(self.spectralMatrixSize);
            for iK=1:size(self.h_pm,2)
                self.h_pm(:,iK) = h(:,self.iK2unique(iK));
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
            newRequiredPropertyNames = {'PFpmInv','QGpmInv','PFpm','QGpm','Ppm','Qpm','QGwg','h_pm'};
        end

        function propertyAnnotations = propertyAnnotationsForGeometry()
            propertyAnnotations = WVGeometryDoublyPeriodicStratified.propertyAnnotationsForGeometry();

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