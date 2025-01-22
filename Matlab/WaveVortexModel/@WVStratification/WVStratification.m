classdef WVStratification < WVAnnotatedClass
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    properties (GetAccess=public, SetAccess=protected)
        % eta_true operation needs rhoFunction
        rhoFunction, N2Function%, dLnN2Function = []
        rho0, rho_nm, N2%, dLnN2
        z, j
        z_int
    end

    properties (Dependent)
        % Z, J % No! these have to be implemented at the transform level
        % because you have to know the full geometry
        Nz, Nj
    end

    methods (Abstract)
        u_z = diffZF(self,u,n);
        w_z = diffZG(self,w,n);
        vm = verticalModes(self);
    end

    properties (Abstract)
        FinvMatrix
        GinvMatrix

        FMatrix
        GMatrix
    end

    methods
        function value = get.Nz(self)
            value=length(self.z);
        end

        function value = get.Nj(self)
            value=length(self.j);
        end

        function flag = isDensityInValidRange(self)
            % checks if the density field is a valid adiabatic re-arrangement of the base state
            %
            % This is probably best re-defined as a dynamical variable.
            %
            % - Topic: Stratification — Validation
            % - Declaration: flag = isDensityInValidRange()
            % - Returns flag: a boolean
            flag = ~(any(self.rho_total(:) < min(self.rho_nm)) | any(self.rho_total(:) > max(self.rho_nm)));
        end

        function effectiveVerticalGridResolution = effectiveVerticalGridResolution(self)
            %returns the effective vertical grid resolution in meters
            %
            % The effective grid resolution is the highest fully resolved
            % wavelength in the model. This value takes into account
            % anti-aliasing, and is thus appropriate for setting damping
            % operators.
            %
            % - Topic: Stratification
            % - Declaration: flag = effectiveVerticalGridResolution(other)
            % - Returns effectiveVerticalGridResolution: double
            arguments
                self WVTransform
            end
            effectiveVerticalGridResolution = pi/max(max(abs(self.l(:)),abs(self.k(:))));
        end

        function throwErrorIfDensityViolation(self,options)
            % checks if the proposed coefficients are a valid adiabatic re-arrangement of the base state
            %
            % Given some proposed new set of values for A0, Ap, Am, will
            % the fluid state violate our density condition? If yes, then
            % throw an error and tell the user about it.
            arguments
                self WVGeostrophicMethods
                options.Ap double = 0
                options.Am double = 0
                options.A0 double = 0
                options.additionalErrorInfo = sprintf(' ');
            end
            % Trying to be extra careful here, so we include the wave part
            % of the flow because we do not really know how everything will
            % add together in the end.
            rho_total = reshape(self.rho_nm,1,1,[]) + (self.rho0/self.g) * shiftdim(self.N2,-2) .* self.transformToSpatialDomainWithG(A0=self.NA0.*options.A0,Apm=self.NAp.*options.Ap + self.NAm.*options.Am);
            densityViolation = any(rho_total(:) < min(self.rho_nm)) | any(rho_total(:) > max(self.rho_nm));
            if densityViolation == 1
                errorString = sprintf('The no-motion density minus rho0 spans from %.3f kg/m^{3} at the surface to %.3f kg/m^{3} at the bottom. Any adiabatic re-arrangement of the fluid requires the density anomaly stay within this range. ',self.rho_nm(end)-self.rho0,self.rho_nm(1)-self.rho0);
                minString = sprintf('\tminimum density: %.3f kg/m^{3}\n',min(rho_total(:))-self.rho0);
                maxString = sprintf('\tmaximum density: %.3f kg/m^{3}\n',max(rho_total(:))-self.rho0);
                errorStruct.message = [errorString,options.additionalErrorInfo,minString,maxString];
                errorStruct.identifier = 'WVTransform:DensityBoundsViolation';
                error(errorStruct);
            end
        end
    end

    methods (Access=protected)
        function self = WVStratification(Lz,Nz,options)
            % By assumption, the modes j=0:(Nj-1) --- because for now, we
            % require completeness.
            % 
            % The required variables for a unique specification of the
            % stratification are (Lz,Nz,N2,Nj,latitude,rho0). However, (Nj,
            % latitude, rho0) all have defaults.
            %
            % To initialize:
            % 1) Pass (Lz,Nz,N2) as a minimum set, and defaults will be
            % chosen (or rho instead of N2).
            % 2) You can also pass everything, and override the defaults,
            % but the assumption is that z is the quadrature grid.
            % 3) verticalModes is strictly optional, and will be created
            % for you if you do not pass it.
            %
            % At the end of initialization, the verticalModes class is
            % initialized with the quadrature point
            arguments
                Lz (1,1) double {mustBePositive}
                Nz (1,1) double {mustBePositive}
                options.z (:,1) double {mustBeNonempty} % quadrature points!
                options.Nj (1,1) double {mustBePositive}
                options.rho function_handle = @isempty
                options.N2 function_handle = @isempty
                options.dLnN2 function_handle = @isempty
                options.latitude (1,1) double = 33
                options.rho0 (1,1) double {mustBePositive} = 1025
                options.verticalModes = []
            end

            % For rigid lid:
            % - There is one barotropic mode that appears in F
            % - There are nModes-1 *internal modes* for G and F.
            % - We compute the nModes+1 internal mode for F, to make it
            % complete.
            % This is nModes+1 grid points necessary to make this happen.
            % This should make sense because there are nModes-1 internal
            % modes, but the boundaries.
            if isfield(options,'z')
                self.z=options.z;
            else
                self.z = WVStratification.quadraturePointsForStratifiedFlow(Lz,Nz,rho=options.rho,N2=options.N2,latitude=options.latitude);
            end
            self.Nz = length(z);

            if isfield(options,'Nj') && isfield(options,'j')
                if options.Nj ~= length(options.j)
                    error('You specified both Nj and j, but they are not compatible!');
                end
            end

            if isfield(options,'Nj')
                if options.Nj > self.Nz-1
                    error('The number of modes must be no greater than Nz-1');
                end
                self.Nj = options.Nj;
            else
                self.Nj = self.Nz-1;
            end

            if self.Nj>1
                self.j = (0:(self.Nj-1))';
            else
                self.j=1;
            end

            nModes = Nz-1;
            if ~isempty(options.verticalModes)
                self.verticalModes = options.verticalModes;
                self.N2Function = options.N2;
                self.rhoFunction = options.rho;
                self.rho_nm = self.rhoFunction(z);
                self.N2 = self.N2Function(z);
            elseif ~isequal(options.N2,@isempty)
                self.verticalModes = InternalModesWKBSpectral(N2=options.N2,zIn=[-Lz 0],zOut=z,latitude=options.latitude,rho0=options.rho0,nModes=nModes,nEVP=max(256,floor(2.1*Nz)));
                self.N2 = options.N2(z);
                self.N2Function = options.N2;
                self.rhoFunction = self.verticalModes.rho_function;
                self.rho_nm = self.rhoFunction(z);
            elseif ~isequal(options.rho,@isempty)
                self.verticalModes = InternalModesWKBSpectral(rho=options.rho,zIn=[-Lz 0],zOut=z,latitude=options.latitude,rho0=options.rho0,nModes=nModes,nEVP=max(256,floor(2.1*Nz)));
                self.N2 = self.verticalModes.N2;
                self.N2Function = self.verticalModes.N2_function;
                self.rhoFunction = options.rho;
                self.rho_nm = self.rhoFunction(z);
            else
                error('You must specify either rho or N2.');
            end
            self.verticalModes.normalization = Normalization.kConstant;
            self.verticalModes.upperBoundary = UpperBoundary.rigidLid;

            if isequal(options.dLnN2,@isempty)
                self.dLnN2 = self.verticalModes.rho_zz./self.verticalModes.rho_z;
            else
                self.dLnN2 = options.dLnN2(z);
                self.dLnN2Function = options.dLnN2;
            end

        end

        function initializeStratifiedFlow(wvt)
            % After initializing the WVTransform, this method can be called
            % and the WVStratification will register.
            arguments
                wvt WVTransform
            end
            wvt.addPropertyAnnotations(WVStratification.propertyAnnotationsForStratifiedFlow);
            wvt.addOperation(EtaTrueOperation());
            wvt.addOperation(APVOperation());
        end

        function [P,Q,PFinv,PF,QGinv,QG,h,w] = verticalProjectionOperatorsForGeostrophicModes(self,Nj)
            % Now go compute the appropriate number of modes at the
            % quadrature points.
            self.verticalModes.normalization = Normalization.geostrophic;
            self.verticalModes.upperBoundary = UpperBoundary.rigidLid;
            [Finv,Ginv,h] = self.verticalModes.ModesAtFrequency(0);
            [P,Q,PFinv,PF,QGinv,QG,h,w] = WVStratification.verticalProjectionOperatorsWithRigidLid(Finv,Ginv,h,Nj,self.verticalModes.Lz);
        end

        function [P,Q,PFinv,PF,QGinv,QG,h] = verticalProjectionOperatorsForIGWModes(self,k,Nj)
            % Now go compute the appropriate number of modes at the
            % quadrature points.
            self.verticalModes.normalization = Normalization.kConstant;
            self.verticalModes.upperBoundary = UpperBoundary.rigidLid;
            [Finv,Ginv,h] = self.verticalModes.ModesAtWavenumber(k);
            [P,Q,PFinv,PF,QGinv,QG,h] = WVStratification.verticalProjectionOperatorsWithRigidLid(Finv,Ginv,h,Nj,self.verticalModes.Lz);
        end
    end

    methods (Static)
        z = quadraturePointsForStratifiedFlow(Lz,Nz,options);
        [P,Q,PFinv,PF,QGinv,QG,h,w] = verticalProjectionOperatorsWithRigidLid(Finv,Ginv,h,Nj,Lz);
        [P,Q,PFinv,PF,QGinv,QG,h] = verticalProjectionOperatorsWithFreeSurface(Finv,Ginv,h,Nj,Lz);

        function writeStratificationToFile(self,ncfile,matFilePath)
            % write the WVStratification to NetCDF and Matlab sidecar file.
            %
            % The NetCDF file must be already initialized and it is assumed
            % that any existing Matlab file at the path is safe to
            % overwrite. This method is designed to be used by the
            % WVTransform classes.
            %
            % % For proper error checking and to write the file
            % independently of the WVTransform classes, use the static
            % method,
            %   `WVStratificationHydrostatic.writeToFile`
            %
            %
            % - Declaration: writeStratificationToFile(ncfile,matFilePath)
            % - Parameter ncfile: a valid NetCDFFile instance
            % - Parameter matFilePath: path to an appropriate location to create a new matlab sidecar file, if needed
            arguments
                self WVStratification {mustBeNonempty}
                ncfile NetCDFFile {mustBeNonempty}
                matFilePath char
            end



            dims = union(dimensions,{'z','j'});
            for iDim=1:length(dims)
                dimAnnotation = dimensionAnnotationNameMap(dims{iDim});
                dimAnnotation.attributes('units') = dimAnnotation.units;
                dimAnnotation.attributes('long_name') = dimAnnotation.description;
                ncfile.addDimension(dimAnnotation.name,self.(dimAnnotation.name),attributes=dimAnnotation.attributes);
            end

            propertyAnnotation = WVStratification.propertyAnnotationsForStratification;
            propertyAnnotationNameMap = configureDictionary("string","WVPropertyAnnotation");
            for i=1:length(propertyAnnotation)
                propertyAnnotationNameMap(propertyAnnotation(i).name) = propertyAnnotation(i);
            end

            requiredVariables = {'rho0','rho_nm','N2'};
            for iVar=1:length(requiredVariables)
                varAnnotation = propertyAnnotationNameMap(requiredVariables{iVar});
                varAnnotation.attributes('units') = varAnnotation.units;
                varAnnotation.attributes('long_name') = varAnnotation.description;
                ncfile.addVariable(varAnnotation.name,varAnnotation.dimensions,self.(varAnnotation.name),isComplex=varAnnotation.isComplex,attributes=varAnnotation.attributes);
            end
        end

        function dimensionAnnotationNameMap = dimensionAnnotationNameMapForStratification(self)
            dimensionAnnotation = WVStratification.dimensionAnnotationsForStratification;
            dimensionAnnotationNameMap = configureDictionary("string","WVDimensionAnnotation");
            for i=1:length(dimensionAnnotation)
                dimensionAnnotationNameMap(dimensionAnnotation(i).name) = dimensionAnnotation(i);
            end
        end


    end
    methods (Static, Hidden=true)
        % All the metadata has to be defined at the class level---so static
        % methods. Thus we have,
        % -dimensionAnnotationsForStratification
        % -propertyAnnotationsForStratification
        % -methodAnnotationsForStratification
        % -requiredDimensionsForStratification
        % -requiredVariablesForStratification
        % function vars = requiredVariablesForStratification()
        %     vars = {'rho0','rho_nm','N2'};
        % end
        % 
        % function dims = requiredDimensionsForStratification()
        %     dims = {'z','j'};
        % end

        function dimensions = dimensionAnnotationsForStratification()
            % return array of WVDimensionAnnotation to annotate the
            % dimensions
            %
            % This function returns annotations for all dimensions of the
            % WVStratification class.
            %
            % - Topic: Developer
            % - Declaration: dimensionAnnotations = WVStratification.dimensionAnnotationsForStratifiedFlow()
            % - Returns dimensionAnnotations: array of WVDimensionAnnotation instances
            dimensions = PMDimensionAnnotation.empty(0,0);

            dimensions(end+1) = PMDimensionAnnotation('z', 'm', 'z coordinate');
            dimensions(end).attributes('standard_name') = 'height_above_mean_sea_level';
            dimensions(end).attributes('positive') = 'up';
            dimensions(end).attributes('axis') = 'Z';

            dimensions(end+1) = PMDimensionAnnotation('j', 'mode number', 'vertical mode number');
        end
        function propertyAnnotations = propertyAnnotationsForStratification()
            % return array of WVPropertyAnnotation initialized by default
            %
            % This function returns annotations for all properties of the
            % WVStratification class.
            %
            % - Topic: Developer
            % - Declaration: propertyAnnotations = WVStratification.propertyAnnotationsForStratifiedFlow()
            % - Returns propertyAnnotations: array of WVPropertyAnnotation instances
            propertyAnnotations = PMPropertyAnnotation.empty(0,0);
            propertyAnnotations(end+1) = PMPropertyAnnotation('verticalModes',{},'', 'instance of the InternalModes class');
            propertyAnnotations(end+1) = PMPropertyAnnotation('rho_nm',{'z'},'kg m^{-3}', '$$\rho_\textrm{nm}(z)$$, no-motion density');
            propertyAnnotations(end+1) = PMPropertyAnnotation('N2',{'z'},'rad^2 s^{-2}', '$$N^2(z)$$, squared buoyancy frequency of the no-motion density, $$N^2\equiv - \frac{g}{\rho_0} \frac{\partial \rho_\textrm{nm}}{\partial z}$$');
            propertyAnnotations(end+1) = PMPropertyAnnotation('dLnN2',{'z'},'', '$$\frac{\partial \ln N^2}{\partial z}$$, vertical variation of the log of the squared buoyancy frequency');
            propertyAnnotations(end+1) = PMPropertyAnnotation('FinvMatrix',{'z','j'},'', 'transformation matrix $$F_g^{-1}$$');
            propertyAnnotations(end+1) = PMPropertyAnnotation('FMatrix',{'j','z'},'', 'transformation matrix $$F_g$$');
            propertyAnnotations(end+1) = PMPropertyAnnotation('GinvMatrix',{'z','j'},'', 'transformation matrix $$G_g^{-1}$$');
            propertyAnnotations(end+1) = PMPropertyAnnotation('GMatrix',{'j','z'},'', 'transformation matrix $$G_g$$');
            propertyAnnotations(end+1) = PMPropertyAnnotation('z_int',{'z'},'', 'Quadrature weights for the vertical grid');
            propertyAnnotations(end+1) = PMPropertyAnnotation('rho0',{},'kg m^{-3}', 'density of $$\rho_\textrm{nm}$$ at the surface (z=0)', detailedDescription='- topic: Domain Attributes');
            propertyAnnotations(end).attributes('standard_name') = 'sea_surface_density';
        end

        % function methodAnnotations = methodAnnotationsForStratification()
        %     % return array of WVAnnotations to annotate the methods
        %     %
        %     % This function returns annotations for all methods of the
        %     % WVStratification class.
        %     %
        %     % - Topic: Developer
        %     % - Declaration: methodAnnotations = WVStratification.methodAnnotationsForStratifiedFlow()
        %     % - Returns methodAnnotations: array of WVAnnotations instances
        %     methodAnnotations = WVAnnotation.empty(0,0);
        % 
        %     methodAnnotations(end+1) = WVAnnotation('diffZF', 'differentiates a variable of (x,y,z) by projecting onto the F-modes, differentiating, and transforming back to (x,y,z)');
        %     methodAnnotations(end+1) = WVAnnotation('diffZG', 'differentiates a variable of (x,y,z) by projecting onto the G-modes, differentiating, and transforming back to (x,y,z)');
        % end
    end

end