classdef WVModel < handle & WVModelAdapativeTimeStepMethods & WVModelFixedTimeStepMethods & WVModelAdapativeTimeStepCellMethods
    % The WVModel is responsible for time-stepping (integrating) the ocean state forward in time, as represented by a WVTransform.
    %
    % Assuming you have already initialized a WVTransform, e.g.,
    % ```matlab
    % wvt = WVTransformConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], N0,latitude=latitude);
    % ```
    % and maybe set some initial conditions, you can then initialize the
    % model,
    % ```matlab
    % model = WVModel(wvt)
    % ```
    % 
    % By default the model only takes a linear time-step. To specify a
    % nonlinear flux on initialization, for example,
    %```matlab
    % model = WVModel(wvt);
    %```
    %
    % You can also initialize a model from existing output,
    % ```matlab
    % model = WVModel.modelFromFile('SomeFile.nc');
    %```
    % 
    % Advanced usage
    % --------------
    %
    % 1. Create 1 or more NetCDFFiles
    % 2. Add 1 or more WVModelOutputGroups to each file
    % 3. Add 1 or more WVObservingSystems to each output group.
    %
    % - Topic: Initialization
    % - Topic: Model Properties
    % - Topic: Integration
    % - Topic: Particles
    % - Topic: Tracer
    % - Topic: Writing to NetCDF files

    properties (GetAccess=public,SetAccess=protected)
        % The WVTransform instance the represents the ocean state.
        % - Topic: Model Properties
        % Set on initialization only, the WVTransform in the model
        % performs all computations necessary to return information about
        % the ocean state at a given time.
        wvt

        fluxedObservingSystems = WVObservingSystem.empty(0,0)
        nFluxComponents = 0
        nFluxComputations uint64 = 0
        indicesForFluxedSystem

        % Indicates whether or not the model is using linear or nonlinear dynamics.
        % - Topic: Model Properties
        % In practice, this is simply checking whether the nonlinearFlux
        % property is nil.
        isDynamicsLinear

        eulerianObservingSystem

        integrationCallback
    end

    properties (Dependent)
        % Current model time (seconds)
        % - Topic: Model Properties
        % Current time of the ocean state, particle positions, and tracer.
        t % (1,1) double

        outputFiles
    end

    properties (Access=private)
        outputFileNameMap = configureDictionary("string","WVModelOutputFile")
    end

    methods (Static)
        model = modelFromFile(path,options)

        function name = defaultOutputGroupName()
            name = "wave-vortex";
        end
    end

    methods

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = WVModel(wvt,options)
            % Initialize a model from a WVTransform instance
            %
            % - Topic: Initialization
            % - Declaration: WVModel(wvt,options)
            % - Parameter wvt: a WaveVortexTranform instance
            %
            % 
            arguments
                wvt WVTransform {mustBeNonempty}
                options.shouldUseLinearDynamics = false
            end

            self.wvt = wvt;
            self.isDynamicsLinear = options.shouldUseLinearDynamics;
            if ~self.isDynamicsLinear
                self.addFluxedCoefficients(WVCoefficients(self));
                if self.wvt.hasClosure == false
                    warning('The nonlinear flux has no damping and may not be stable.');
                end
            end
            self.eulerianObservingSystem = WVEulerianFields(self,fieldNames=intersect({'Ap','Am','A0'},self.wvt.variableNames));
        end

        function value = get.t(self)
            value = self.wvt.t;
        end
        
        function ncfile = ncfile(self)
            ncfile = NetCDFFile.empty(0,0);
            outputFiles_ = self.outputFiles;
            for iFile = 1:length(outputFiles_)
                if ~isempty(outputFiles_(iFile).ncfile)
                    ncfile(end+1) = outputFiles_(iFile).ncfile;
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Output groups
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function outputFiles = get.outputFiles(self)
            outputFiles = [self.outputFileNameMap(self.outputFileNameMap.keys)];
        end
        function names = outputFileNames(self)
            % retrieve the names of all output group names
            %
            % - Topic: Utility function — Metadata
            arguments (Input)
                self WVModel {mustBeNonempty}
            end
            arguments (Output)
                names string
            end
            names = self.outputFileNameMap.keys;
        end

        function val = outputFileWithName(self,name)
            % retrieve a WVModelOutputGroup by name
            arguments (Input)
                self WVModel {mustBeNonempty}
                name char {mustBeNonempty}
            end
            arguments (Output)
                val WVModelOutputGroup
            end
            val = self.outputFileNameMap(name);
        end

        function addOutputFile(self,outputFile)
            arguments
                self WVModel {mustBeNonempty}
                outputFile WVModelOutputFile
            end
            self.outputFileNameMap(outputFile.filename) = outputFile;
        end

        function outputFile = addNewOutputFile(self,path,options)
            arguments (Input)
                self WVModel {mustBeNonempty}
                path {mustBeText}
                options.shouldOverwriteExisting logical = false
            end
            arguments (Output)
                outputFile WVModelOutputFile
            end
            outputFile = WVModelOutputFile(self,path,shouldOverwriteExisting=options.shouldOverwriteExisting);
            self.addOutputFile(outputFile);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Fluxed observing systems
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function wvCoeff = wvCoefficientFluxedObservingSystem(self)
            wvCoeff = [];
            if ~isempty(self.fluxedObservingSystems) && isa(self.fluxedObservingSystems(1),'WVCoefficients')
                wvCoeff = self.fluxedObservingSystems(1);
            end
        end

        function addFluxedCoefficients(self,anObservingSystem)
            arguments
                self WVModel {mustBeNonempty}
                anObservingSystem (1,1) WVCoefficients
            end
            % prepend, so that its always first
            if isempty(self.fluxedObservingSystems)
                self.fluxedObservingSystems = anObservingSystem;
            else
                if ~isa(self.fluxedObservingSystems(1),'WVCoefficients')
                    idx = find(size(self.fluxedObservingSystems) > 1, 1);
                    if isempty(idx)
                        idx = 1;
                    end
                    self.fluxedObservingSystems = cat(idx, anObservingSystem, self.fluxedObservingSystems);
                end
            end
            self.recomputeIndicesForFluxedSystems();
        end 

        function addFluxedObservingSystem(self,anObservingSystem)
            arguments
                self WVModel {mustBeNonempty}
                anObservingSystem WVObservingSystem
            end
            for iObs=1:length(anObservingSystem)
                alreadyInArray = any(cellfun(@(x) isequal(x, anObservingSystem(iObs)), num2cell(self.fluxedObservingSystems)));
                if ~alreadyInArray
                    self.fluxedObservingSystems(end+1) = anObservingSystem(iObs);
                end
            end
            self.recomputeIndicesForFluxedSystems();
        end

        function removeFluxedObservingSystem(self,anObservingSystem)
            arguments
                self WVModel {mustBeNonempty}
                anObservingSystem WVObservingSystem
            end
            for iObs=1:length(anObservingSystem)
                self.fluxedObservingSystems = setdiff(self.fluxedObservingSystems,anObservingSystem(iObs),'stable');
            end
            self.recomputeIndicesForFluxedSystems();
        end

        function anObservingSystem = fluxedObservingSystemWithName(self,name)
            arguments (Input)
                self WVModel {mustBeNonempty}
                name {mustBeText}
            end
            arguments (Output)
                anObservingSystem WVObservingSystem
            end
            anObservingSystem = [];
            for i = 1:length(self.fluxedObservingSystems)
                if strcmp(self.fluxedObservingSystems(i).name,name)
                    anObservingSystem = self.fluxedObservingSystems(i);
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Convenience functions for adding/removing Eulerian variables
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function addNetCDFOutputVariables(self,variables)
            % Add variables to list of variables to be written to the NetCDF variable during the model run.
            %
            % - Topic: Writing to NetCDF files
            % - Declaration: addNetCDFOutputVariables(variables)
            % - Parameter variables: strings of variable names.
            %
            % Pass strings of WVTransform state variables of the
            % same name. This must be called before using any of the
            % integrate methods.
            %
            % ```matlab
            % model.addNetCDFOutputVariables('A0','u','v');
            % ```
            
            arguments
                self WVModel
            end
            arguments (Repeating)
                variables char
            end
            self.eulerianObservingSystem.addNetCDFOutputVariables(variables{:});
        end

        function setNetCDFOutputVariables(self,variables)
            % Set list of variables to be written to the NetCDF variable during the model run.
            %
            % - Topic: Writing to NetCDF files
            % - Declaration: setNetCDFOutputVariables(variables)
            % - Parameter variables: strings of variable names.
            %
            % Pass strings of WVTransform state variables of the
            % same name. This must be called before using any of the
            % integrate methods.
            %
            % ```matlab
            % model.setNetCDFOutputVariables('A0','u','v');
            % ```
            arguments
                self WVModel
            end
            arguments (Repeating)
                variables char
            end
            self.eulerianObservingSystem.setNetCDFOutputVariables(variables{:});
        end

        function removeNetCDFOutputVariables(self,variables)
            % Remove variables from the list of variables to be written to the NetCDF variable during the model run.
            %
            % - Topic: Writing to NetCDF files
            % - Declaration: removeNetCDFOutputVariables(variables)
            % - Parameter variables: strings of variable names.
            %
            % Pass strings of WVTransform state variables of the
            % same name. This must be called before using any of the
            % integrate methods.
            %
            % ```matlab
            % model.removeNetCDFOutputVariables('A0','u','v');
            % ```
            arguments
                self WVModel
            end
            arguments (Repeating)
                variables char
            end
            self.eulerianObservingSystem.removeNetCDFOutputVariables(variables{:});
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Convenience functions for floats and drifters and tracer
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function addParticles(self,name,isXYOnly,x,y,z,trackedFieldNames,options)
            % Add particles to be advected by the flow.
            %
            % - Topic: Particles
            % - Declaration: addParticles(name,fluxOp,x,y,z,trackedFieldNames,options)
            % - Parameter name: a unique name to call the particles
            % - Parameter fluxOp: a WVParticleFluxOperation, used to determine how the flow advects the particles
            % - Parameter x: x-coordinate location of the particles
            % - Parameter y: y-coordinate location of the particles
            % - Parameter z: z-coordinate location of the particles
            % - Parameter trackedFields: strings of variable names
            % - Parameter advectionInterpolation: (optional) interpolation method used for particle advection. "linear" (default), "spline", "exact"
            % - Parameter trackedVarInterpolation: (optional) interpolation method used for tracked field. "linear" (default), "spline", "exact"
            % - Parameter absToleranceXY: (adapative) absolute tolerance in meters for particle advection in (x,y). 1e-1 (default)
            % - Parameter absToleranceZ: (adapative) absolute tolerance  in meters for particle advection in (z). 1e-2 (default)
            arguments
                self WVModel {mustBeNonempty}
                name char {mustBeNonempty}
                isXYOnly logical {mustBeNonempty}
                x (1,:) double
                y (1,:) double
                z (1,:) double
            end
            arguments (Repeating)
                trackedFieldNames char
            end
            arguments
                options.advectionInterpolation char {mustBeMember(options.advectionInterpolation,["linear","spline","exact","finufft"])} = "linear"
                options.trackedVarInterpolation char {mustBeMember(options.trackedVarInterpolation,["linear","spline","exact","finufft"])} = "spline"
                options.outputGroupName = "wave-vortex"
                options.absToleranceXY = 1e-1; % 100 km * 10^{-6}
                options.absToleranceZ = 1e-2;
            end

            observingSystem = WVLagrangianParticles(self,name=name,isXYOnly=isXYOnly,x=x,y=y,z=z,trackedFieldNames=trackedFieldNames,advectionInterpolation=options.advectionInterpolation,trackedVarInterpolation=options.trackedVarInterpolation);
            if isscalar(self.outputFiles) && isscalar(self.outputFiles(1).outputGroups)
                self.outputFiles(1).outputGroups(1).addObservingSystem(observingSystem);
            elseif isempty(self.outputFiles)
                self.addFluxedObservingSystem(observingSystem);
            else
                error('There is more than one output file associated with this model. You must manually choose which file to add particles to.');
            end
        end

        function [x,y,z,trackedFields] = particlePositions(self,name)
            % Positions and values of tracked fields of particles at the current model time.
            %
            % - Topic: Particles
            % - Declaration: [x,y,z,trackedFields] = particlePositions(name)
            % - Parameter name: name of the particles
            [x,y,z,trackedFields] = self.fluxedObservingSystemWithName(name).particlePositions;
        end
        
        function setFloatPositions(self,x,y,z,trackedFields,options)
            % Set positions of float-like particles to be advected by the model.
            %
            % - Topic: Particles
            % - Declaration: setFloatPositions(self,x,y,z,trackedFields,options)
            % - Parameter x: x-coordinate location of the particles
            % - Parameter y: y-coordinate location of the particles
            % - Parameter z: z-coordinate location of the particles
            % - Parameter trackedFields: strings of variable names
            % - Parameter advectionInterpolation: (optional) interpolation method used for particle advection. "linear" (default), "spline", "exact"
            % - Parameter trackedVarInterpolation: (optional) interpolation method used for tracked field. "linear" (default), "spline", "exact"
            % - Parameter absToleranceXY: (adapative) absolute tolerance in meters for particle advection in (x,y). 1e-1 (default)
            % - Parameter absToleranceZ: (adapative) absolute tolerance  in meters for particle advection in (z). 1e-2 (default)
            %
            % Pass the initial positions of particles to be advected by all
            % three components of the velocity field, (u,v,w).
            %
            % Particles move between grid (collocation) points and thus
            % their location must be interpolated. By default the
            % advectionInterpolation is set to "linear" interpolation. For
            % many flows this will have sufficient accuracy and allow you
            % to place float at nearly every grid point without slowing
            % down the model integration. However, if high accuracy is
            % required, you may want to use cubic "spline" interpolation or
            % even "exact" at the expense of computational speed.
            %
            % You can track the value of any known WVVariableAnnotation along the
            % particle's flow path, e.g., relative vorticity. These values
            % must also be interpolated using one of the known
            % interpolation methods.
            %
            % ```matlab
            % nTrajectories = 101;
            % xFloat = Lx/2*ones(1,nTrajectories);
            % yFloat = Ly/2*ones(1,nTrajectories);
            % zFloat = linspace(-Lz,0,nTrajectories);
            %
            % model.setFloatPositions(xFloat,yFloat,zFloat,'rho_total');
            % ```
            %
            % If a NetCDF file is set for output, the particle positions
            % and tracked fields will automatically be written to file
            % during integration. If you are not writing to file you can
            % retrieve the current positions and values of the tracked
            % fields by calling -floatPositions.
            arguments
                self WVModel {mustBeNonempty}
                x (1,:) double
                y (1,:) double
                z (1,:) double = []
            end
            arguments (Repeating)
                trackedFields char
            end
            arguments
                options.advectionInterpolation char {mustBeMember(options.advectionInterpolation,["linear","spline","exact","finufft"])} = "linear"
                options.trackedVarInterpolation char {mustBeMember(options.trackedVarInterpolation,["linear","spline","exact","finufft"])} = "linear"
                options.absToleranceXY = 1e-1;
                options.absToleranceZ = 1e-2;
            end
            optionCell = namedargs2cell(options);
            self.addParticles('float',false,x,y,z,trackedFields{:},optionCell{:});
        end

        function [x,y,z,tracked] = floatPositions(self)
            % Returns the positions of the floats at the current time as well as the value of the fields being tracked.
            %
            % - Topic: Particles
            % - Declaration: [x,y,z,tracked] = floatPositions()
            %
            % The tracked variable is a structure, with fields named for
            % each of the requested fields being tracked.
            %
            % In the following example, float positions are set along with
            % one tracked field 
            % ```matlab
            % model.setFloatPositions(xFloat,yFloat,zFloat,'rho_total');
            %
            % % Set up the integrator
            % nT = model.setupIntegrator(timeStepConstraint="oscillatory", outputInterval=period/10,finalTime=3*period);
            %
            % % write the float trajectories to memory
            % xFloatT = zeros(nT,nTrajectories);
            % yFloatT = zeros(nT,nTrajectories);
            % zFloatT = zeros(nT,nTrajectories);
            % rhoFloatT = zeros(nT,nTrajectories);
            % t = zeros(nT,1);
            %
            % [xFloatT(1,:),yFloatT(1,:),zFloatT(1,:),tracked] = model.floatPositions;
            % rhoFloatT(1,:) = tracked.rho_total;
            % ```
            %
            [x,y,z,tracked] = self.particlePositions('float');
        end

        function setDrifterPositions(self,x,y,z,trackedFields,options)
            % Set positions of drifter-like particles to be advected.
            % - Topic: Particles
            arguments
                self WVModel {mustBeNonempty}
                x (1,:) double
                y (1,:) double
                z (1,:) double = []
            end
            arguments (Repeating)
                trackedFields char
            end
            arguments
                options.advectionInterpolation char {mustBeMember(options.advectionInterpolation,["linear","spline","exact","finufft"])} = "linear"
                options.trackedVarInterpolation char {mustBeMember(options.trackedVarInterpolation,["linear","spline","exact","finufft"])} = "linear"
                options.absToleranceXY = 1e-1;
            end
            optionCell = namedargs2cell(options);
            self.addParticles('drifter',true,x,y,z,trackedFields{:},optionCell{:});
        end

        function [x,y,z,tracked] = drifterPositions(self)
            % Current positions of the drifter particles
            % - Topic: Particles
            [x,y,z,tracked] = self.particlePositions('drifter');
        end

        function addTracer(self,phi,name)
            arguments
                self WVModel
                phi double
                name {mustBeNonempty,mustBeText}
            end
            % Add a scalar field tracer to be advected by the flow
            % - Topic: Tracer
            isXYOnly= (length(self.wvt.spatialDimensionNames) == 2);
            observingSystem = WVTracer(self,name=name,phi=phi,isXYOnly=isXYOnly);
            if isscalar(self.outputFiles) && isscalar(self.outputFiles(1).outputGroups)
                self.outputFiles(1).outputGroups(1).addObservingSystem(observingSystem);
            elseif isempty(self.outputFiles)
                self.addFluxedObservingSystem(observingSystem);
            else
                error('There is more than one output file associated with this model. You must manually choose which file to add particles to.');
            end
        end

        function phi = tracer(self,name)
            % Scalar field of the requested tracer at the current model time.
            % - Topic: Tracer
            phi = self.fluxedObservingSystemWithName(name).phi;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Integration loop
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function setupIntegrator(self,options,fixedTimeStepOptions,adaptiveTimeStepOptions)
            % Customize the time-stepping
            %
            % By default the model will use adaptive time stepping with a
            % reasonable choice of values. However, you may find it
            % necessary to customize the time stepping behavior.
            %
            % When setting up the integrator you must choice between
            % "adaptive" and "fixed" integrator types. Depending on which
            % type you choose, you will have different options available.
            %
            % The "fixed" time-step integrator used a cfl condition based
            % on the advective velocity, but you can change this to use the
            % highest oscillatory frequency. Alternatively, you can simply
            % set deltaT yourself.
            %
            % The "adaptive" time-step integator uses absolute and relative
            % error tolerances. It is worth reading Matlab's documentation
            % on RelTol and AbsTol as part of odeset to understand what
            % these mean. By default, the adaptive time stepping uses a
            % a relative error tolerance of 1e-3 for everything. However,
            % the absolute error tolerance is less straightforward.
            %
            % The absolute tolerance has a meaningful scale with units, and
            % thus must be chosen differently for particle positions (x,y)
            % than for geostrophic coefficients (A0). 
            %
            % - Topic: Integration
            % - Declaration: setupIntegrator(self,options)
            % - Parameter integratorType: (optional) integrator type. "adaptive"(default), "fixed"
            % - Parameter deltaT: (fixed) time step
            % - Parameter cfl: (fixed) cfl condition
            % - Parameter timeStepConstraint: (fixed) constraint to fix the time step. "advective" (default) ,"oscillatory","min"
            % - Parameter integrator: (adapative) function handle of integrator. @ode78 (default)
            % - Parameter absTolerance: (adapative) absolute tolerance for sqrt(energy). 1e-6 (default)
            % - Parameter relTolerance: (adapative) relative tolerance for sqrt(energy). 1e-3 (default)
            % - Parameter shouldShowIntegrationStats: (adapative) whether to show integration output 0 or 1 (default)
            arguments
                self WVModel {mustBeNonempty}

                options.integratorType char {mustBeMember(options.integratorType,["fixed","adaptive","adaptive-cell"])} = "adaptive"

                fixedTimeStepOptions.deltaT (1,1) double {mustBePositive}
                fixedTimeStepOptions.cfl (1,1) double
                fixedTimeStepOptions.timeStepConstraint char {mustBeMember(fixedTimeStepOptions.timeStepConstraint,["advective","oscillatory","min"])} = "min"

                adaptiveTimeStepOptions.integrator = @ode78
                adaptiveTimeStepOptions.absTolerance = 1e-6
                adaptiveTimeStepOptions.relTolerance = 1e-3;
                adaptiveTimeStepOptions.shouldShowIntegrationStats double {mustBeMember(adaptiveTimeStepOptions.shouldShowIntegrationStats,[0 1])} = 0
            end

            if self.isDynamicsLinear == false
                self.wvCoefficientFluxedObservingSystem.absTolerance = adaptiveTimeStepOptions.absTolerance;
            end

            % self.resetFixedTimeStepIntegrator();
            self.resetAdapativeTimeStepIntegrator();

            self.integratorType = options.integratorType;
            if strcmp(self.integratorType,"adaptive")
                adaptiveTimeStepOptions = rmfield(adaptiveTimeStepOptions,"absTolerance");
                optionArgs = namedargs2cell(adaptiveTimeStepOptions);
                self.setupAdaptiveTimeStepIntegrator(optionArgs{:});
            elseif strcmp(self.integratorType,"adaptive-cell")
                adaptiveTimeStepOptions = rmfield(adaptiveTimeStepOptions,"absTolerance");
                adaptiveTimeStepOptions = rmfield(adaptiveTimeStepOptions,"integrator");
                optionArgs = namedargs2cell(adaptiveTimeStepOptions);
                self.setupAdaptiveTimeStepCellIntegrator(optionArgs{:});
            else
                optionArgs = namedargs2cell(fixedTimeStepOptions);
                self.setupFixedTimeStepIntegrator(optionArgs{:});
            end 

            self.didSetupIntegrator = true;
        end
        
        
        function integrateToTime(self,finalTime,options)
            % Time step the model forward to the requested time.
            % - Topic: Integration
            arguments
                self WVModel {mustBeNonempty}
                finalTime (1,:) double
                options.shouldShowIntegrationDiagnostics logical = true
                options.shouldAllowBackwardsIntegration = false
                options.callback
            end
            if finalTime <= self.t && ~options.shouldAllowBackwardsIntegration
                fprintf('Reqested integration to time %d, but the model is currently at time t=%d.\n',round(finalTime),round(self.t));
                return;
            end
            if ~self.didSetupIntegrator
                self.setupIntegrator();
            end

            % if self.nFluxComponents == 0
            %     if self.eulerianObservingSystem.nTimeSeriesVariables == 0
            %         error("Nothing to do! There are no variables being integrated and no dynamical fields being output.");
            %     else
            %         warning('There no variables being integrated, and the variables that are being written can be recovered instantly from the initial conditions.');
            %     end
            % end

            self.shouldShowIntegrationDiagnostics = options.shouldShowIntegrationDiagnostics;
            if isfield(options,'callback')
                self.integrationCallback = options.callback;
            end
                  
            % arrayfun( @(outputFile) outputFile.initializeOutputFile(), self.outputFiles);
            % arrayfun( @(outputFile) outputFile.writeTimeStepToOutputFile(self.t), self.outputFiles);

            self.wvt.restoreForcingAmplitudes();
            

            if self.nFluxComponents == 0
                self.pseudoIntegrateToTime(finalTime);
            elseif strcmp(self.integratorType,"adaptive")
                self.integrateToTimeWithAdaptiveTimeStep(finalTime)
            elseif strcmp(self.integratorType,"adaptive-cell")
                self.integrateToTimeWithAdaptiveTimeStepCell(finalTime)
            else
                self.integrateToTimeWithFixedTimeStep(finalTime);
            end

            self.recordNetCDFFileHistory();            
        end


        function recomputeIndicesForFluxedSystems(self)
            self.nFluxComponents = 0;
            for i = 1:length(self.fluxedObservingSystems)
                self.nFluxComponents = self.fluxedObservingSystems(i).nFluxComponents + self.nFluxComponents;
            end

            self.indicesForFluxedSystem = cell(length(self.fluxedObservingSystems),1);
            nFinal = 0;
            for i = 1:length(self.fluxedObservingSystems)
                nInitial = nFinal + 1;
                nFinal = nFinal + self.fluxedObservingSystems(i).nFluxComponents;
                
                self.indicesForFluxedSystem{i} = reshape(nInitial:nFinal,[],1);
            end
        end

        function Y0 = absErrorToleranceCellArray(self)
            Y0 = cell(self.nFluxComponents,1);
            for i = 1:length(self.fluxedObservingSystems)
                Y0(self.indicesForFluxedSystem{i}) = self.fluxedObservingSystems(i).absErrorTolerance();
            end
        end

        function Y0 = initialConditionsCellArray(self)
            Y0 = cell(self.nFluxComponents,1);
            for i = 1:length(self.fluxedObservingSystems)
                Y0(self.indicesForFluxedSystem{i}) = self.fluxedObservingSystems(i).initialConditions();
            end
        end

        function F = fluxAtTimeCellArray(self,t,y0)
            self.nFluxComputations = self.nFluxComputations + 1;
            F = cell(self.nFluxComponents,1);
            for i = 1:length(self.fluxedObservingSystems)
                F(self.indicesForFluxedSystem{i}) = self.fluxedObservingSystems(i).fluxAtTime(t,y0(self.indicesForFluxedSystem{i}));
            end
        end

        function updateIntegratorValuesFromCellArray(self,t,y0)
            % We must set the time here. If we are integrating the
            % wave-vortex coefficients, then this is benign because it will
            % immediately get repeated momentarily. But we we are not
            % integrating, and are running linearly, then we need the
            % fields to update.
            self.wvt.t = t;
            for i = 1:length(self.fluxedObservingSystems)
                self.fluxedObservingSystems(i).updateIntegratorValues(t,y0(self.indicesForFluxedSystem{i}));
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % NetCDF Output
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function summarize(self)
            if isempty(self.fluxedObservingSystems)
                fprintf('The model has no systems to integrate.\n');
            else
                if isscalar(self.fluxedObservingSystems)
                    fprintf("The model is integrating 1 system.\n");
                else
                    fprintf("The model is integrating " + length(self.fluxedObservingSystems) + " systems.\n");
                end

                for i = 1:length(self.fluxedObservingSystems)
                    fprintf("\t" + string(i) + ": " + self.fluxedObservingSystems(i).description + "\n");
                end
            end
            fprintf('\n');
            if isempty(length(self.outputFiles))
                fprintf('The model has no output files.\n');
            else
                if isscalar(self.outputFiles)
                    fprintf("The model will write to 1 output file.\n");
                else
                    fprintf("The model will write to " + length(self.outputFiles) + " output files.\n");
                end
                for iFile=1:length(self.outputFiles)
                    if isscalar(length(self.outputFiles(iFile).outputGroups))
                        fprintf(string(iFile) + ": " + self.outputFiles(iFile).filename + " with " + length(self.outputFiles(iFile).outputGroups) + " output group\n");
                    else
                        fprintf(string(iFile) + ": " + self.outputFiles(iFile).filename + " with " + length(self.outputFiles(iFile).outputGroups) + " output groups\n");
                    end
                    for iGroup=1:length(self.outputFiles(iFile).outputGroups)
                        nObsSystems = length(self.outputFiles(iFile).outputGroups(iGroup).observingSystems);
                        fprintf("\t" + string(iGroup) + ": " + self.outputFiles(iFile).outputGroups(iGroup).description + ", writing " + string(nObsSystems) + " observing systems\n");
                        for iOs = 1:nObsSystems
                            fprintf("\t\t" + string(iOs) + ": " + self.outputFiles(iFile).outputGroups(iGroup).observingSystems(iOs).description + "\n");
                        end
                    end
                end
            end
        end

        function outputFile = createNetCDFFileForModelOutput(self,path,options)
            % Create a NetCDF file for model output
            % - Topic: Writing to NetCDF files
            arguments
                self WVModel {mustBeNonempty}
                path char {mustBeNonempty}
                options.outputInterval (1,1) double {mustBePositive}
                options.shouldOverwriteExisting logical = false
            end

            outputFile = self.addNewOutputFile(path,shouldOverwriteExisting=options.shouldOverwriteExisting);
            outputGroup = outputFile.addNewEvenlySpacedOutputGroup(self.defaultOutputGroupName,initialTime=self.t,outputInterval=options.outputInterval);

            outputGroup.addObservingSystem(self.eulerianObservingSystem);
            for i = 1:length(self.fluxedObservingSystems)
                outputGroup.addObservingSystem(self.fluxedObservingSystems(i));
            end 

            %% Now what happens?
            % Still need to set the default group? No, its set by name
            % Do we need to trigger a write? Or let the first time step do
            % that?
        end

        function recordNetCDFFileHistory(self,options)
            arguments
                self WVModel {mustBeNonempty}
                options.didBlowUp {mustBeNumeric} = 0
            end

            arrayfun( @(outputFile) outputFile.recordNetCDFFileHistory(didBlowUp=options.didBlowUp), self.outputFiles);
        end

    end


    properties %(Access = protected)
        didSetupIntegrator=false

        integrationStartWallTime
        integrationStartModelTime
        integrationLastInformWallTime       % wall clock, to keep track of the expected integration time
        integrationLastInformModelTime
        integrationInformTime = 10
        nFluxComputationsAtLastInform uint64 = 0

        integratorType      % Array integrator
        finalIntegrationTime % set only during an integration
        shouldShowIntegrationDiagnostics = true

        % Initial model time (seconds)
        % - Topic: Model Properties
        % The time of the WVTransform when the model was
        % initialized. This also corresponds to the first time in the
        % NetCDF output file.
        initialTime (1,1) double = 0 
    end
    

    methods %(Access=protected)

        function flag = didBlowUp(self)
            if ( any(isnan(self.wvt.Ap)|isnan(self.wvt.Am)|isnan(self.wvt.A0)) )
                flag = 1;
                fprintf('Blowup detected. Aborting.');
            else
                flag = 0;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Diagnostics (start/during/finish) integration
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function showIntegrationStartDiagnostics(self,finalTime)
            if self.shouldShowIntegrationDiagnostics  == false
                return;
            end
            self.nFluxComputations = 0;
            self.integrationStartWallTime = datetime('now');
            self.integrationStartModelTime = self.wvt.t;
            fprintf('Starting numerical simulation on %s.\n', datetime(self.integrationStartWallTime,TimeZone='local',Format='d-MMM-y HH:mm:ss Z'));
            fprintf('\tStarting at model time t=%.2f inertial periods and integrating to t=%.2f inertial periods.\n',self.t/self.wvt.inertialPeriod,finalTime/self.wvt.inertialPeriod);
            self.integrationLastInformWallTime = datetime('now');
            self.integrationLastInformModelTime = self.wvt.t;
        end

        function showIntegrationTimeDiagnostics(self,finalTime)
            if self.shouldShowIntegrationDiagnostics == false && isempty(self.integrationCallback)  == 0
                return;
            end
            deltaWallTime = datetime('now')-self.integrationLastInformWallTime;
            if ( seconds(deltaWallTime) > self.integrationInformTime)
                if self.shouldShowIntegrationDiagnostics
                    wallTimePerModelTime = deltaWallTime / (self.wvt.t - self.integrationLastInformModelTime);
                    wallTimeRemaining = wallTimePerModelTime*(finalTime - self.wvt.t);
                    deltaT = (self.wvt.t-self.integrationLastInformModelTime)/( self.nFluxComputations - self.nFluxComputationsAtLastInform);
                    fprintf('\tmodel time t=%.2f inertial periods. Estimated time to reach %.2f inertial periods is %s (%s). Δ≅%.2fs\n', self.t/self.wvt.inertialPeriod, finalTime/self.wvt.inertialPeriod, wallTimeRemaining, datetime(datetime('now')+wallTimeRemaining,TimeZone='local',Format='d-MMM-y HH:mm:ss Z'),deltaT) ;
                    self.wvt.summarizeEnergyContent();

                    self.integrationLastInformWallTime = datetime('now');
                    self.integrationLastInformModelTime = self.wvt.t;
                    self.nFluxComputationsAtLastInform = self.nFluxComputations;
                end
                if ~isempty(self.integrationCallback)
                    self.integrationCallback(self);
                end
            end
        end

        function showIntegrationFinishDiagnostics(self)
            self.integrationCallback = [];
            if self.shouldShowIntegrationDiagnostics  == false
                return;
            end
            integrationTotalTime = datetime('now')-self.integrationStartWallTime;
            deltaT = (self.wvt.t-self.integrationStartModelTime)/self.nFluxComputations ;
            fprintf('Finished after time %s. Δ≅%.2fs\n', integrationTotalTime,deltaT);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Write to file
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function pseudoIntegrateToTime(self,finalTime)
            % Time step the model forward linearly
            arguments
                self WVModel {mustBeNonempty}
                finalTime (1,1) double
            end

            % The function call here is stupid, because it is not obvious
            % that callign outputTimesForIntegrationPeriod actually has the
            % side-effect of setting up the run
            integratorTimes = self.outputTimesForIntegrationPeriod(self.t,finalTime);
            arrayfun( @(outputFile) outputFile.writeTimeStepToOutputFile(self.t), self.outputFiles);

            self.finalIntegrationTime = finalTime;
            for iTime=1:length(integratorTimes)
                self.wvt.t = integratorTimes(iTime);
                self.writeTimeStepToNetCDFFile(self.wvt.t);
            end
            self.finalIntegrationTime = [];
        end

        function integratorTimes = outputTimesForIntegrationPeriod(self,initialTime,finalTime)
            % This will be called exactly once before an integration
            % begins.
            arguments (Input)
                self WVModel
                initialTime (1,1) double
                finalTime (1,1) double
            end
            arguments (Output)
                integratorTimes (:,1) double
            end
            integratorTimes = [];
            outputFiles_ = self.outputFiles;
            for iFile = 1:length(outputFiles_)
                integratorTimes = cat(1,integratorTimes,outputFiles_(iFile).outputTimesForIntegrationPeriod(initialTime,finalTime));
            end
            integratorTimes = sort(uniquetol(integratorTimes));
            if isempty(integratorTimes) || integratorTimes(1) ~= self.t
                integratorTimes = cat(1,self.t,integratorTimes);
            end
            if integratorTimes(end) ~= finalTime
                integratorTimes = cat(1,integratorTimes,finalTime);
            end
        end

        function writeTimeStepToNetCDFFile(self,t)
            outputFiles_ = self.outputFiles;
            for iFile = 1:length(outputFiles_)
                outputFiles_(iFile).writeTimeStepToOutputFile(t);
            end
        end

        function closeNetCDFFile(self)
            arrayfun( @(outputFile) outputFile.closeNetCDFFile(), self.outputFiles);
        end
    end


end