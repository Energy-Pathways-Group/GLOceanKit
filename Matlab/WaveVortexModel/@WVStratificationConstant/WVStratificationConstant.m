classdef WVStratificationConstant < WVStratification
    properties (GetAccess=public, SetAccess=protected)
        N0
        h_pm
        h_0
    end

    properties (Dependent)
        isHydrostatic
        FinvMatrix
        GinvMatrix
        FMatrix
        GMatrix
    end

    methods
        function self = WVStratificationConstant(Lz, Nz, options)
            % create matrices for hydrostatics in variable stratification
            %
            % To initialize:
            % 1) Pass (Lz,Nz,N0) as a minimum set, and defaults will be
            % chosen.
            % 2) Pass (Lz,Nz,N0,Nj,latitude,rho0) to fully, uniquely
            % specify the transforms
            %
            % - Topic: Initialization
            % - Declaration: wvt = WVTransformHydrostatic(Lxyz, Nxyz, options)
            % - Parameter Lxyz: length of the domain (in meters) in the three coordinate directions, e.g. [Lx Ly Lz]
            % - Parameter Nxyz: number of grid points in the three coordinate directions, e.g. [Nx Ny Nz]
            % - Parameter rho:  (optional) function_handle specifying the density as a function of depth on the domain [-Lz 0]
            % - Parameter stratification:  (optional) function_handle specifying the stratification as a function of depth on the domain [-Lz 0]
            % - Parameter latitude: (optional) latitude of the domain (default is 33 degrees north)
            % - Parameter rho0: (optional) density at the surface z=0 (default is 1025 kg/m^3)
            % - Returns wvt: a new WVTransformHydrostatic instance
            arguments
                Lz (1,1) double {mustBePositive}
                Nz (1,1) double {mustBePositive}
                options.N0 (1,1) double {mustBePositive} = 5.2e-3
                options.Nj (1,1) double {mustBePositive}
                options.latitude (1,1) double = 33
                options.rho0 (1,1) double {mustBePositive} = 1025
            end
            if isfield(options,'Nj')
                if options.Nj > self.Nz-1
                    error('The number of modes must be no greater than Nz-1');
                end
                Nj = options.Nj;
            else
                Nj = Nz-1;
            end
            dz = Lz/(options.Nz-1);
            self.z = dz*(0:(options.Nz-1))' - Lz; % Cosine basis for DCT-I and DST-I
            self.z_int = dz*ones(Nz,1);
            self.z_int(1) = dz/2; self.z_int(end) = dz/2;
            self.j = (0:(Nj-1))';
            self.N0 = options.N0;
            self.rho0 = options.rho0;
            self.rhoFunction = @(z) -(N0*N0*rho0/9.81)*z + rho0;
            self.N2Function = @(z) N0*N0*ones(size(z));
            self.rho_nm = self.rhoFunction(self.z);
            self.N2 = self.N2Function(self.z);            
        end

        function vm = verticalModes(self)
            vm = InternalModesConstantStratification(N0=self.N0, rho0=self.rho0, zIn=[self.z(1) 0], zOut=self.z, latitude=self.latitude);
        end

        function writeStratifiedFlowToFile(self,ncfile,matFilePath)
            % write the WVStratificationHydrostatic to NetCDF and Matlab sidecar file.
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
            % - Declaration: writeStratifiedFlowToFile(ncfile,matFilePath)
            % - Parameter ncfile: a valid NetCDFFile instance
            % - Parameter matFilePath: path to an appropriate location to create a new matlab sidecar file, if needed
            arguments
                self WVStratificationConstant {mustBeNonempty}
                ncfile NetCDFFile {mustBeNonempty}
                matFilePath char
            end
            % This will add the dimensions to the root of the file
            writeStratifiedFlowToFile@WVStratifiedFlow(self,ncfile,matFilePath);

            % To keep things tidy, lets put the transform pieces in a separate group
            group = ncfile.addGroup("WVStratificationConstant");

            propertyAnnotation = WVStratificationConstant.defaultPropertyAnnotations;
            propertyAnnotationNameMap = configureDictionary("string","WVPropertyAnnotation");
            for i=1:length(propertyAnnotation)
                propertyAnnotationNameMap(propertyAnnotation(i).name) = propertyAnnotation(i);
            end

            requiredVariables = {'N0'};
            for iVar=1:length(requiredVariables)
                varAnnotation = propertyAnnotationNameMap(requiredVariables{iVar});
                varAnnotation.attributes('units') = varAnnotation.units;
                varAnnotation.attributes('long_name') = varAnnotation.description;
                group.addVariable(varAnnotation.name,varAnnotation.dimensions,self.(varAnnotation.name),isComplex=varAnnotation.isComplex,attributes=varAnnotation.attributes);
            end
        end

        function initializeStratifiedFlow(wvt)
            % After initializing the WVTransform, this method can be called
            % and the WVStratifiedFlow will register.
            arguments
                wvt WVTransform
            end
            wvt.addPropertyAnnotations(WVStratificationConstant.propertyAnnotationsForStratifiedFlow);
            wvt.addOperation(EtaTrueOperation());
            wvt.addOperation(APVOperation());
        end
    end

    methods (Access=protected)
        function vars = requiredVariablesForStratification(self)
            vars = {'N0'};
        end
        function dims = requiredDimensionsForStratification(self)
            dims = {'z','j'};
        end
    end

    methods (Static)
        % All the metadata has to be defined at the class level---so static
        % methods. Thus we have,
        % -dimensionAnnotationsForStratification
        % -propertyAnnotationsForStratification
        % -methodAnnotationsForStratification
        % -requiredDimensionsForStratification
        % -requiredVariablesForStratification
        function propertyAnnotations = propertyAnnotationsForStratification()
            % return array of WVPropertyAnnotation initialized by default
            %
            % This function returns annotations for all properties of the
            % WVStratificationConstant class (as well as its
            % superclass).
            %
            % - Topic: Internal
            % - Declaration: propertyAnnotations = WVStratificationConstant.propertyAnnotationsForStratification()
            % - Returns propertyAnnotations: array of WVPropertyAnnotation instances
            propertyAnnotations = WVStratification.propertyAnnotationsForStratification();
            propertyAnnotations(end+1) = WVPropertyAnnotation('N0',{},'rad s^{-1}', 'buoyancy frequency of the no-motion density');
        end

        function [bool, errorString] = canInitializeDirectlyFromFile(ncfile)
            bool = false;
            errorString = "";
            requiredDimensions = WVStratificationConstant.requiredDimensions();
            if ~all(ncfile.hasVariableWithName(requiredDimensions{:}))
                errorString = "The NetCDF file is missing a required dimension.";
                return;
            end

            if ~ncfile.hasGroupWithName("WVStratificationConstant")
                errorString = "The NetCDF is missing the group WVStratificationHydrostatic.";
                return;
            end

            group = ncfile.groupWithName("WVStratificationConstant");
            requiredVariables = WVStratificationConstant.requiredVariables();
            if ~all(group.hasVariableWithName(requiredVariables{:}))
                errorString = "The NetCDF group WVStratificationConstant is missing a required variable.";
                return;
            end

            bool = true;
        end

        function var = requiredVariablesFromFile(ncfile)
            requiredDimensions = WVStratificationConstant.requiredDimensions();
            for iDim = 1:length(requiredDimensions)
                name = requiredDimensions{iDim};
                var.(name) = ncfile.readVariables(name);
            end
            requiredVariables = WVStratificationConstant.requiredVariables();
            for iDim = 1:length(requiredVariables)
                name = requiredVariables{iDim};
                var.(name) = ncfile.readVariables(name);
            end
        end

        % 1) instance method to write to NetCDF group
        % 2) static method declaring all required variables
        % 3) init methods takes group as argument, reads required variables
        % 4) static method can read netcdf file
        function flow = stratifiedFlowFromFile(group,wvt)
            arguments
                group NetCDFGroup {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end
            flow = WVSpectralVanishingViscosity(wvt,nu_xy=group.attributes('nu_xy'),nu_z=group.attributes('nu_z') );
        end
        
    end
end