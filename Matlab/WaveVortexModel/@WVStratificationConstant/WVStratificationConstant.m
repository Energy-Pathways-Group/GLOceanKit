classdef WVStratificationConstant < WVStratification
    properties (GetAccess=public, SetAccess=protected)
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
            % 1) Pass (Lz,Nz,N2) as a minimum set, and defaults will be
            % chosen (or rho instead of N2).
            % 2) Pass (Lz,Nz,N2,Nj,latitude,rho0) to fully, uniquely
            % specify the transforms
            % 3) 
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
            group = ncfile.addGroup("WVStratificationHydrostatic");

            propertyAnnotation = WVStratificationConstant.defaultPropertyAnnotations;
            propertyAnnotationNameMap = configureDictionary("string","WVPropertyAnnotation");
            for i=1:length(propertyAnnotation)
                propertyAnnotationNameMap(propertyAnnotation(i).name) = propertyAnnotation(i);
            end

            requiredVariables = {'dLnN2','PF0inv','QG0inv','PF0','QG0','P0','Q0','h','z_int'};
            for iVar=1:length(requiredVariables)
                varAnnotation = propertyAnnotationNameMap(requiredVariables{iVar});
                varAnnotation.attributes('units') = varAnnotation.units;
                varAnnotation.attributes('long_name') = varAnnotation.description;
                group.addVariable(varAnnotation.name,varAnnotation.dimensions,self.(varAnnotation.name),isComplex=varAnnotation.isComplex,attributes=varAnnotation.attributes);
            end
            
            rhoFunction = self.rhoFunction;
            N2Function = self.N2Function;
            dLnN2Function = self.dLnN2Function;
            date_created = ncfile.attributes('date_created');
            save(matFilePath,'rhoFunction','N2Function','dLnN2Function','date_created');
            fprintf('In addition to the NetCDF file, a .mat sidecar file was created at the same path.\n');
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

    methods (Static)
        function propertyAnnotations = propertyAnnotationsForStratifiedFlow()
            % return array of WVPropertyAnnotation initialized by default
            %
            % This function returns annotations for all properties of the
            % WVStratificationHydrostatic class (as well as its
            % superclass).
            %
            % - Topic: Internal
            % - Declaration: propertyAnnotations = WVStratificationHydrostatic.propertyAnnotationsForStratifiedFlow()
            % - Returns propertyAnnotations: array of WVPropertyAnnotation instances
            propertyAnnotations = WVStratifiedFlow.propertyAnnotationsForStratifiedFlow();
            propertyAnnotations(end+1) = WVPropertyAnnotation('PF0inv',{'z','j'},'','Preconditioned F-mode inverse transformation');
            propertyAnnotations(end+1) = WVPropertyAnnotation('QG0inv',{'z','j'},'','Preconditioned G-mode inverse transformation');
            propertyAnnotations(end+1) = WVPropertyAnnotation('PF0',{'j','z'},'','Preconditioned F-mode forward transformation');
            propertyAnnotations(end+1) = WVPropertyAnnotation('QG0',{'j','z'},'','Preconditioned G-mode forward transformation');
            propertyAnnotations(end+1) = WVPropertyAnnotation('P0',{'j'},'','Preconditioner for F, size(P)=[1 Nj]. F*u = uhat, (PF)*u = P*uhat, so ubar==P*uhat');
            propertyAnnotations(end+1) = WVPropertyAnnotation('Q0',{'j'},'','Preconditioner for G, size(Q)=[1 Nj]. G*eta = etahat, (QG)*eta = Q*etahat, so etabar==Q*etahat. ');
            propertyAnnotations(end+1) = WVPropertyAnnotation('h',{'j'},'m', 'equivalent depth of each mode', detailedDescription='- topic: Domain Attributes — Stratification');
            propertyAnnotations(end+1) = WVPropertyAnnotation('z_int',{'z'},'', 'Quadrature weights for the vertical grid');
        end

        function vars = requiredVariables()
            vars = {'dLnN2','PF0inv','QG0inv','PF0','QG0','P0','Q0','h','z_int'};
        end

        function dims = requiredDimensions()
            dims = {'z','j'};
        end

        function [ncfile,matFilePath] = writeToFile(stratifiedFlow,path,options)
            % Output the WVStratificationHydrostatic to file.
            %
            % Writes the WVStratificationHydrostatic instance to file, with enough information to
            % re-initialize.
            %
            % - Topic: Write to file
            % - Declaration: [ncfile,matFilePath] = WVStratificationHydrostatic.writeToFile(stratifiedFlow,path,options)
            % - Parameter stratifiedFlow: instance of WVStratificationHydrostatic
            % - Parameter path: path to write the file.
            % - Parameter shouldOverwriteExisting: (optional) boolean indicating whether or not to overwrite an existing file at the path. Default 0.
            arguments (Input)
                stratifiedFlow WVStratificationConstant {mustBeNonempty}
                path char {mustBeNonempty}
            end
            arguments (Input)
                options.shouldOverwriteExisting double {mustBeMember(options.shouldOverwriteExisting,[0 1])} = 0
            end
            arguments (Output)
                ncfile NetCDFFile
                matFilePath char
            end

            matFilePath = WVStratifiedFlow.matlabSidecarPathForNetCDFPath(path);

            if options.shouldOverwriteExisting == 1
                if isfile(path)
                    delete(path);
                end
                if isfile(matFilePath)
                    delete(matFilePath);
                end
            else
                if isfile(path) || isfile(matFilePath)
                    error('A file already exists with that name.')
                end
            end
            ncfile = NetCDFFile(path);
            stratifiedFlow.writeStratifiedFlowToFile(ncfile,matFilePath);
        end

        function [bool, errorString] = canInitializeDirectlyFromFile(ncfile)
            bool = false;
            errorString = "";
            matFilePath = WVStratifiedFlow.matlabSidecarPathForNetCDFPath(ncfile.path);
            if ~isfile(matFilePath)
                errorString = "The .mat sidecar file is missing, which is necessary for direct initialization.";
                return;
            end
            matFile = load(matFilePath);
            if isa(matFile.N2Function,'function_handle') == false && isa(matFile.rhoFunction,'function_handle') == false
                errorString = "The .mat sidecar file does not contain either an N2 or a rho function_handle object.";
                return;
            end
            requiredDimensions = WVStratificationConstant.requiredDimensions();
            if ~all(ncfile.hasVariableWithName(requiredDimensions{:}))
                errorString = "The NetCDF file is missing a required dimension.";
                return;
            end

            if ~ncfile.hasGroupWithName("WVStratificationHydrostatic")
                errorString = "The NetCDF is missing the group WVStratificationHydrostatic.";
                return;
            end

            group = ncfile.groupWithName("WVStratificationHydrostatic");
            requiredVariables = WVStratificationConstant.requiredVariables();
            if ~all(group.hasVariableWithName(requiredVariables{:}))
                errorString = "The NetCDF group WVStratificationHydrostatic is missing a required variable.";
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

            matFilePath = WVStratifiedFlow.matlabSidecarPathForNetCDFPath(ncfile.path);
            matFile = load(matFilePath);
            if isa(matFile.N2Function,'function_handle') == true
                var.N2 = matFile.N2Function;
            else
                var.N2 = @isempty;
            end
            if isa(matFile.rhoFunction,'function_handle') == true
                var.rho = matFile.rhoFunction;
            else
                var.rho = @isempty;
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