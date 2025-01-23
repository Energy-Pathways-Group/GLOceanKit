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

        function ncfile = writeStratificationToFile(self,path,options)
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
            arguments (Input)
                self WVStratificationConstant {mustBeNonempty}
                path char {mustBeNonempty}
                options.shouldOverwriteExisting double {mustBeMember(options.shouldOverwriteExisting,[0 1])} = 0
            end
            arguments (Output)
                ncfile NetCDFFile
            end
            options.dims = WVStratificationConstant.requiredDimensionsForStratification();
            options.properties = WVStratificationConstant.requiredPropertiesForStratification();
            options.dimAnnotations = WVStratificationConstant.dimensionAnnotationsForStratification();
            options.propAnnotations = WVStratificationConstant.propertyAnnotationsForStratification();
            namedOptions = namedargs2cell(options);
            ncfile = PMAnnotatedClass.writeToPath(self,path,namedOptions{:});
        end

        function initializeStratifiedFlow(wvt)
            % After initializing the WVTransform, this method can be called
            % and the WVStratifiedFlow will register.
            arguments
                wvt WVTransform
            end
            % wvt.addPropertyAnnotations(WVStratificationConstant.propertyAnnotationsForStratifiedFlow);
            wvt.addOperation(EtaTrueOperation());
            wvt.addOperation(APVOperation());
        end
    end

    methods (Access=protected)

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
            propertyAnnotations(end+1) = PMPropertyAnnotation('N0',{},'rad s^{-1}', 'buoyancy frequency of the no-motion density');
        end

        function vars = requiredPropertiesForStratification()
            vars = {'N0','latitude','rho0'};
        end
        function dims = requiredDimensionsForStratification()
            dims = {'z','j'};
        end

        % 1) instance method to write to NetCDF group
        % 2) static method declaring all required variables
        % 3) init methods takes group as argument, reads required variables
        % 4) static method can read netcdf file
        function stratification = stratificationFromFile(path)
            arguments (Input)
                path char {mustBeNonempty}
            end
            arguments (Output)
                stratification WVStratificationConstant {mustBeNonempty}
            end
            ncfile = NetCDFFile(path);
            stratification = WVStratificationConstant.stratificationFromGroup(ncfile);
        end

        function stratification = stratificationFromGroup(group)
            arguments (Input)
                group NetCDFGroup {mustBeNonempty}
            end
            arguments (Output)
                stratification WVStratificationConstant {mustBeNonempty}
            end
            requiredVariables = WVStratificationConstant.requiredPropertiesForStratification;
            requiredDimensions = WVStratificationConstant.requiredDimensionsForStratification;
            [canInit, errorString] = canInitializeDirectlyFromGroup(group,requiredDimensions,requiredVariables);
            if ~canInit
                error(errorString);
            end

            requiredVariables = union(requiredVariables,requiredDimensions);
            vars = PMAnnotatedClass.variablesFromGroup(group,requiredVariables);
            
            Nz = length(vars.z);
            Lz = vars.z(end) - vars.z(1);
            options.Nj = length(vars.j);
            options.latitude = vars.latitude;
            options.rho0 = vars.rho0;
            options.N0 = vars.N0;
            optionCell = namedargs2cell(options);
            stratification = WVStratificationConstant(Lz,Nz,optionCell{:});
        end
        
    end
end