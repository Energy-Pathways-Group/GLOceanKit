classdef WVStratificationHydrostatic < WVStratification
    properties (Hidden=true) %(GetAccess=public, SetAccess=protected) %(Access=private)
        % Transformation matrices
        PF0inv, QG0inv % size(PFinv,PGinv)=[Nz x Nj]
        PF0, QG0 % size(PF,PG)=[Nj x Nz]
        h % [Nj 1]

        P0 % Preconditioner for F, size(P)=[Nj 1]. F*u = uhat, (PF)*u = P*uhat, so ubar==P*uhat
        Q0 % Preconditioner for G, size(Q)=[Nj 1]. G*eta = etahat, (QG)*eta = Q*etahat, so etabar==Q*etahat.

        zInterp
        PFinvInterp, QGinvInterp
    end

    properties (GetAccess=private, SetAccess=private)
        Nj
    end

    properties (Dependent)
        isHydrostatic
        FinvMatrix
        GinvMatrix
        FMatrix
        GMatrix
    end

    methods
        function self = WVStratificationHydrostatic(Lz, Nz, options)
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
                options.z (:,1) double {mustBeNonempty} % quadrature points!
                options.Nj (1,1) double {mustBePositive}
                options.rho function_handle = @isempty
                options.N2 function_handle = @isempty
                options.dLnN2 function_handle = @isempty
                options.latitude (1,1) double = 33
                options.rho0 (1,1) double {mustBePositive} = 1025
                options.verticalModes = []
                options.ncfile NetCDFFile
                options.matFilePath
            end
            canInitializeDirectly = false;
            if isfield(options,'ncfile')
                requiredDimensions = {'z','j'};
                requiredVariables = {'dLnN2','PF0inv','QG0inv','PF0','QG0','P0','Q0','h','z_int'};
                if ncfile.hasVariableWithName(requiredDimensions{:}) && ncfile.hasGroupWithName("WVStratificationHydrostatic")
                    group = ncfile.groupWithName("WVStratificationHydrostatic");
                    if group.hasVariableWithName(requiredVariables{:})
                        canInitializeDirectly = true;
                    end
                end
            end

            if canInitializeDirectly == true

            else
                [self.P0,self.Q0,self.PF0inv,self.PF0,self.QG0inv,self.QG0,self.h,self.z_int] = self.verticalProjectionOperatorsForGeostrophicModes(self.Nj);
            end

            % First we need to initialize the WVStratifiedFlow.
            if isfield(options,'z')
                z = options.z;
            else
                z = WVStratifiedFlow.quadraturePointsForStratifiedFlow(Lz,Nz,rho=options.rho,N2=options.N2,latitude=options.latitude);
            end
            self@WVStratifiedFlow(Lz,z,rho=options.rho,N2=options.N2,dLnN2=options.dLnN2func,latitude=options.latitude)

            % if all of these things are set initially (presumably read
            % from file), then we can initialize without computing modes.
            canInitializeDirectly = all(isfield(options,{'N2','latitude','rho0','dLnN2','PFinv','QGinv','PF','QG','h','P','Q','z'}));

            if canInitializeDirectly
                fprintf('Initialize the WVTransformHydrostatic directly from matrices.\n');
                self.Nj = size(options.PF,1);
            else
                nModes = Nxyz(3)-1;
                if options.shouldAntialias == 1
                    self.Nj = floor(options.jAliasingFraction*nModes);
                else
                    self.Nj = nModes;
                end
            end

            if canInitializeDirectly
                self.PF0inv = options.PFinv;
                self.QG0inv = options.QGinv;
                self.PF0 = options.PF;
                self.QG0 = options.QG;
                self.h = options.h;
                self.P0 = options.P;
                self.Q0 = options.Q;
            else
                [self.P0,self.Q0,self.PF0inv,self.PF0,self.QG0inv,self.QG0,self.h,self.z_int] = self.verticalProjectionOperatorsForGeostrophicModes(self.Nj);
            end


            % self.offgridModes = WVOffGridTransform(im,self.latitude, self.N2Function,1);

            self.dftBuffer = zeros(self.spatialMatrixSize);
            self.wvBuffer = zeros([self.Nz self.Nkl]);
            [self.dftPrimaryIndex, self.dftConjugateIndex, self.wvConjugateIndex] = self.horizontalModes.indicesFromWVGridToDFTGrid(self.Nz,isHalfComplex=1);
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
                self WVStratificationHydrostatic {mustBeNonempty}
                ncfile NetCDFFile {mustBeNonempty}
                matFilePath char
            end
            % This will add the dimensions to the root of the file
            writeStratifiedFlowToFile@WVStratifiedFlow(self,ncfile,matFilePath);

            % To keep things tidy, lets put the transform pieces in a separate group
            group = ncfile.addGroup("WVStratificationHydrostatic");

            propertyAnnotation = WVStratificationHydrostatic.defaultPropertyAnnotations;
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
            wvt.addPropertyAnnotations(WVStratificationHydrostatic.propertyAnnotationsForStratifiedFlow);
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
                stratifiedFlow WVStratificationHydrostatic {mustBeNonempty}
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
            requiredDimensions = WVStratificationHydrostatic.requiredDimensions();
            if ~all(ncfile.hasVariableWithName(requiredDimensions{:}))
                errorString = "The NetCDF file is missing a required dimension.";
                return;
            end

            if ~ncfile.hasGroupWithName("WVStratificationHydrostatic")
                errorString = "The NetCDF is missing the group WVStratificationHydrostatic.";
                return;
            end

            group = ncfile.groupWithName("WVStratificationHydrostatic");
            requiredVariables = WVStratificationHydrostatic.requiredVariables();
            if ~all(group.hasVariableWithName(requiredVariables{:}))
                errorString = "The NetCDF group WVStratificationHydrostatic is missing a required variable.";
                return;
            end

            bool = true;
        end

        function var = requiredVariablesFromFile(ncfile)
            requiredDimensions = WVStratificationHydrostatic.requiredDimensions();
            for iDim = 1:length(requiredDimensions)
                name = requiredDimensions{iDim};
                var.(name) = ncfile.readVariables(name);
            end
            requiredVariables = WVStratificationHydrostatic.requiredVariables();
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