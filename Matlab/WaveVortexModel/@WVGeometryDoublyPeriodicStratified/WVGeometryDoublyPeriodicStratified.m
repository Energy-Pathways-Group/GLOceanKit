classdef WVGeometryDoublyPeriodicStratified < WVGeometryDoublyPeriodic & WVStratification & WVGeometryCartesianXYZ
    properties (Access=public) %(GetAccess=private, SetAccess=private) %(Access=private)
        dLnN2

        % Transformation matrices
        PF0inv, QG0inv % size(PFinv,PGinv)=[Nz x Nj]
        PF0, QG0 % size(PF,PG)=[Nj x Nz]
        h_0 % [Nj 1]
        h_pm

        P0 % Preconditioner for F, size(P)=[Nj 1]. F*u = uhat, (PF)*u = P*uhat, so ubar==P*uhat
        Q0 % Preconditioner for G, size(Q)=[Nj 1]. G*eta = etahat, (QG)*eta = Q*etahat, so etabar==Q*etahat.
    end

    properties (Dependent)
        FinvMatrix
        GinvMatrix
        FMatrix
        GMatrix
        Lr2
        kPseudoRadial
          % [Nj 1]
    end

    methods
        function self = WVGeometryDoublyPeriodicStratified(Lxyz, Nxyz, geomOptions, stratOptions, directInit)
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
                stratOptions.z (:,1) double {mustBeNonempty} % quadrature points!
                stratOptions.j (:,1) double {mustBeNonempty}
                stratOptions.Nj (1,1) double {mustBePositive}
                stratOptions.rhoFunction function_handle = @isempty
                stratOptions.N2Function function_handle = @isempty
                stratOptions.rho0 (1,1) double {mustBePositive} = 1025
                stratOptions.planetaryRadius (1,1) double = 6.371e6
                stratOptions.rotationRate (1,1) double = 7.2921E-5
                stratOptions.latitude (1,1) double = 33
                stratOptions.g (1,1) double = 9.81

                % ALL of these must be set for direct initialization to
                % avoid actually computing the modes.
                directInit.dLnN2 (:,1) double
                directInit.PF0inv
                directInit.QG0inv
                directInit.PF0
                directInit.QG0
                directInit.h_0 (:,1) double
                directInit.P0 (:,1) double
                directInit.Q0 (:,1) double
                directInit.z_int (:,1) double
            end
            Lz = Lxyz(3);
            Nz = Nxyz(3);
            if ~isfield(stratOptions,'z')
                stratOptions.z = WVStratification.quadraturePointsForStratifiedFlow(Lz,Nz,rho=stratOptions.rhoFunction,N2=stratOptions.N2Function,latitude=stratOptions.latitude,rotationRate=stratOptions.rotationRate);
            end
            nModes = Nz-1;
            if ~isequal(stratOptions.N2Function,@isempty)
                verticalModes = InternalModesWKBSpectral(N2=stratOptions.N2Function,zIn=[-Lz 0],zOut=stratOptions.z,latitude=stratOptions.latitude,rho0=stratOptions.rho0,nModes=nModes,nEVP=max(256,floor(2.1*Nz)),rotationRate=stratOptions.rotationRate,g=stratOptions.g);
                stratOptions.N2Function = stratOptions.N2Function;
                stratOptions.rhoFunction = @(z) verticalModes.rho_function(z);
            elseif ~isequal(stratOptions.rhoFunction,@isempty)
                verticalModes = InternalModesWKBSpectral(rho=stratOptions.rhoFunction,zIn=[-Lz 0],zOut=stratOptions.z,latitude=stratOptions.latitude,rho0=stratOptions.rho0,nModes=nModes,nEVP=max(256,floor(2.1*Nz)),rotationRate=stratOptions.rotationRate,g=stratOptions.g);
                stratOptions.N2Function = @(z) verticalModes.N2_function(z);
                stratOptions.rhoFunction = stratOptions.rhoFunction;
            end
            verticalModes.normalization = Normalization.kConstant;
            verticalModes.upperBoundary = UpperBoundary.rigidLid;

            if geomOptions.shouldAntialias == true && ~isfield(stratOptions,"Nj")
                maxNj = Nxyz(3)-1;
                if maxNj > 3
                    stratOptions.Nj = floor(2*maxNj/3);
                end
            end

            statOptionCell = namedargs2cell(stratOptions);
            self@WVStratification(Lxyz(3),Nxyz(3),statOptionCell{:});

            optionCell = namedargs2cell(geomOptions);
            self@WVGeometryDoublyPeriodic(Lxyz(1:2),Nxyz(1:2),optionCell{:},Nz=Nxyz(3),shouldExcludeNyquist=true,shouldExludeConjugates=true,conjugateDimension=2);

            self.verticalModes = verticalModes;

            allFields = cell2struct([struct2cell(stratOptions);struct2cell(directInit)],[fieldnames(stratOptions);fieldnames(directInit)]);
            canInitializeDirectly = all(isfield(allFields, WVGeometryDoublyPeriodicStratified.namesOfRequiredPropertiesForStratification));

            if canInitializeDirectly == true
                self.dLnN2 = directInit.dLnN2;
                self.PF0inv = directInit.PF0inv;
                self.QG0inv = directInit.QG0inv;
                self.PF0 = directInit.PF0;
                self.QG0 = directInit.QG0;
                self.P0 = directInit.P0;
                self.Q0 = directInit.Q0;
                self.h_0 = directInit.h_0;
                self.z_int = directInit.z_int;
            else
                self.dLnN2 = self.verticalModes.rho_zz./self.verticalModes.rho_z;
                [self.P0,self.Q0,self.PF0inv,self.PF0,self.QG0inv,self.QG0,self.h_0,self.z_int] = self.verticalProjectionOperatorsForGeostrophicModes(self.Nj);
            end
            self.h_pm = self.h_0;
        end

        du = diffZF(self,u,options)
        dw = diffZG(self,w,options)

        function Lr2 = get.Lr2(self)
            Lr2 = self.g*self.h_0/(self.f*self.f);
        end

        function j_max = effectiveJMax(self)
            j_max = max(self.j);
        end

        function kPseudoRadial = get.kPseudoRadial(self)
            jWavenumber = 1./sqrt(self.Lr2);
            jWavenumber(1) = 0; % barotropic mode is a mean?
            [kj,kr] = ndgrid(jWavenumber,self.kRadial);
            Kh = sqrt(kj.^2 + kr.^2);
            allKs = unique(reshape(abs(Kh),[],1),'sorted');
            deltaK = max(diff(allKs));
            kAxis_ = 0:deltaK:(max(allKs)+deltaK/2);
            % This choices of axis spacing ensures that there will be no
            % gaps in the resulting spectrum.
            kPseudoRadial = reshape(kAxis_,[],1);
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
                wvt         WVStratification
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
                wvt         WVStratification
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
                wvt         WVStratification
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
                wvt         WVStratification
            end
            G = wvt.QG0 ./ shiftdim(wvt.Q0,2);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Transformations FROM the spatial domain
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function u = transformToSpatialDomainWithF(self, options)
            arguments
                self WVTransform {mustBeNonempty}
                options.Apm double = 0
                options.A0 double = 0
            end
            u = self.transformToSpatialDomainWithFourier(self.PF0inv*(self.P0 .* (options.Apm + options.A0)));
        end

        function w = transformToSpatialDomainWithG(self, options)
            arguments
                self WVTransform {mustBeNonempty}
                options.Apm double = 0
                options.A0 double = 0
            end
            w = self.transformToSpatialDomainWithFourier(self.QG0inv*(self.Q0 .* (options.Apm + options.A0)));
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Spectra
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        S_f = spectrumWithFgTransform(self,f)
        S_f = spectrumWithGgTransform(self,f)
        S_f = crossSpectrumWithFgTransform(self,phi,gamma)
        S_f = crossSpectrumWithGgTransform(self,phi,gamma)

        [varargout] = transformToPseudoRadialWavenumber(self,energyReservoir,varargin);
        [varargout] = transformToPseudoRadialWavenumberA0(self,varargin);
        [varargout] = transformToPseudoRadialWavenumberApm(self,varargin)  
        [varargout] = transformToOmegaAxis(self,varargin)   

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Transformations TO the spatial domain
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function u_bar = transformFromSpatialDomainWithFio(self, u)
            u_bar = (self.PF0*u)./self.P0;
        end

        function u_bar = transformFromSpatialDomainWithFg(self, u)
            u_bar = (self.PF0*u)./self.P0;
        end

        function w_bar = transformFromSpatialDomainWithGg(self, w)
            w_bar = (self.QG0*w)./self.Q0;
        end

        function w_bar = transformWithG_wg(~, w_bar )
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
    end

    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % CAAnnotatedClass required methods, which enables writeToFile
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = WVGeometryDoublyPeriodicStratified.propertyAnnotationsForGeometry();
        end

        function vars = classRequiredPropertyNames()
            vars = WVGeometryDoublyPeriodicStratified.namesOfRequiredPropertiesForGeometry();
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Stratification specific property annotations and initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function requiredPropertyNames = namesOfRequiredPropertiesForGeometry()
            requiredPropertyNames = WVStratification.namesOfRequiredPropertiesForStratification();
            requiredPropertyNames = union(requiredPropertyNames,WVGeometryDoublyPeriodicStratified.newRequiredPropertyNames());
            requiredPropertyNames = union(requiredPropertyNames,WVGeometryDoublyPeriodic.namesOfRequiredPropertiesForGeometry());
            requiredPropertyNames = setdiff(requiredPropertyNames,WVGeometryDoublyPeriodicStratified.newNonrequiredPropertyNames);
        end

        function newRequiredPropertyNames = newRequiredPropertyNames()
            newRequiredPropertyNames = {'dLnN2','PF0inv','QG0inv','PF0','QG0','P0','Q0','h_0','z_int'};
        end

        function newNonrequiredPropertyNames = newNonrequiredPropertyNames()
            newNonrequiredPropertyNames = {'conjugateDimension','shouldExcludeNyquist','shouldExludeConjugates'};
        end

        function propertyAnnotations = propertyAnnotationsForGeometry()
            propertyAnnotations = WVGeometryDoublyPeriodic.propertyAnnotationsForGeometry();
            propertyAnnotations = cat(2,propertyAnnotations,WVStratification.propertyAnnotationsForStratification());
            propertyAnnotations = cat(2,propertyAnnotations,WVGeometryCartesianXYZ.propertyAnnotationsForGeometry());

            propertyAnnotations(end+1) = CANumericProperty('PF0inv',{'z','j'},'','Preconditioned F-mode inverse transformation');
            propertyAnnotations(end+1) = CANumericProperty('QG0inv',{'z','j'},'','Preconditioned G-mode inverse transformation');
            propertyAnnotations(end+1) = CANumericProperty('PF0',{'j','z'},'','Preconditioned F-mode forward transformation');
            propertyAnnotations(end+1) = CANumericProperty('QG0',{'j','z'},'','Preconditioned G-mode forward transformation');
            propertyAnnotations(end+1) = CANumericProperty('P0',{'j'},'','Preconditioner for F, size(P)=[1 Nj]. F*u = uhat, (PF)*u = P*uhat, so ubar==P*uhat');
            propertyAnnotations(end+1) = CANumericProperty('Q0',{'j'},'','Preconditioner for G, size(Q)=[1 Nj]. G*eta = etahat, (QG)*eta = Q*etahat, so etabar==Q*etahat. ');

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
            [Lxyz(3), Nxyz(3), stratOptions] = WVStratification.requiredPropertiesForStratificationFromGroup(group);
            vars = CAAnnotatedClass.propertyValuesFromGroup(group,WVGeometryDoublyPeriodicStratified.newRequiredPropertyNames);
            newOptions = namedargs2cell(vars);
            options = cat(2,stratOptions,geomOptions,newOptions);
        end

        function geometry = geometryFromFile(path)
            arguments (Input)
                path char {mustBeNonempty}
            end
            arguments (Output)
                geometry WVGeometryDoublyPeriodicStratified {mustBeNonempty}
            end
            ncfile = NetCDFFile(path);
            geometry = WVGeometryDoublyPeriodicStratified.geometryFromGroup(ncfile);
        end

        function geometry = geometryFromGroup(group)
            arguments (Input)
                group NetCDFGroup {mustBeNonempty}
            end
            arguments (Output)
                geometry WVGeometryDoublyPeriodicStratified {mustBeNonempty}
            end
            CAAnnotatedClass.throwErrorIfMissingProperties(group,WVGeometryDoublyPeriodicStratified.namesOfRequiredPropertiesForGeometry);
            [Lxyz, Nxyz, options] = WVGeometryDoublyPeriodicStratified.requiredPropertiesForGeometryFromGroup(group);
            geometry = WVGeometryDoublyPeriodicStratified(Lxyz,Nxyz,options{:});
        end

    end
end