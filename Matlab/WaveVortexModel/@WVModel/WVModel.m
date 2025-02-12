classdef WVModel < handle & WVModelAdapativeTimeStepMethods & WVModelFixedTimeStepMethods
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
        
        % Reference to the NetCDFFile being used for model output
        % - Topic: Writing to NetCDF files
        % Empty indicates no file output.
        ncfile NetCDFFile

        % output index of the current/most recent step.
        % - Topic: Integration
        % If stepsTaken=0, outputIndex=1 means the initial conditions get written at index 1
        incrementsWrittenToFile (1,1) uint64 = 0

        % output index of the current/most recent step.
        % - Topic: Integration
        % If stepsTaken=0, outputIndex=1 means the initial conditions get written at index 1
        timeOfLastIncrementWrittenToFile (1,1) double = -Inf
    end

    properties (Dependent)
        % Indicates whether or not the model is using linear or nonlinear dynamics.
        % - Topic: Model Properties
        % In practice, this is simply checking whether the nonlinearFlux
        % property is nil.
        linearDynamics

        % Current model time (seconds)
        % - Topic: Model Properties
        % Current time of the ocean state, particle positions, and tracer.
        t % (1,1) double

        outputGroups
    end

    properties (Access=private)
        outputGroupNameMap = configureDictionary("string","WVModelOutputGroup")
    end

    methods (Static)
        model = modelFromFile(path,options)

        function name = defaultOutputGroupName()
            name = "wave-vortex";
        end
    end

    methods
        function outputGroups = get.outputGroups(self)
            outputGroups = [self.outputGroupNameMap(self.outputGroupNameMap.keys)];
        end
        function names = outputGroupNames(self)
            % retrieve the names of all output group names
            %
            % - Topic: Utility function — Metadata
            arguments (Input)
                self WVModel {mustBeNonempty}
            end
            arguments (Output)
                names string
            end
            names = self.outputGroupNameMap.keys;
        end

        function val = outputGroupWithName(self,name)
            % retrieve a WVModelOutputGroup by name
            arguments (Input)
                self WVModel {mustBeNonempty}
                name char {mustBeNonempty}
            end
            arguments (Output)
                val WVModelOutputGroup
            end
            val = self.outputGroupNameMap(name);
        end

        function outputGroup = addOutputGroup(self,name,options)
            arguments
                self WVModel {mustBeNonempty}
                name {mustBeText}
                options.outputInterval
            end
            if isfield(options,"outputInterval")
                outputGroup = WVModelOutputGroup(self,name,outputInterval=options.outputInterval);
            else
                outputGroup = WVModelOutputGroup(self,name);
            end
            self.outputGroupNameMap(name) = outputGroup;
        end

        function addNetCDFOutputVariables(self,variables,options)
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
            arguments
                options.outputGroupName
            end
            if ~isfield(options,"outputGroupName")
                options.outputGroupName = self.defaultOutputGroupName;
            end
            outputGroup = self.outputGroupWithName(options.outputGroupName);
            outputGroup.addNetCDFOutputVariables(variables{:});
        end

        function setNetCDFOutputVariables(self,variables,options)
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
            arguments
                options.outputGroupName
            end
            if ~isfield(options,"outputGroupName")
                options.outputGroupName = self.defaultOutputGroupName;
            end
            outputGroup = self.outputGroupWithName(options.outputGroupName);
            outputGroup.setNetCDFOutputVariables(variables{:});
        end

        function removeNetCDFOutputVariables(self,variables,options)
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
            arguments
                options.outputGroupName
            end
            if ~isfield(options,"outputGroupName")
                options.outputGroupName = self.defaultOutputGroupName;
            end
            outputGroup = self.outputGroupWithName(options.outputGroupName);
            outputGroup.removeNetCDFOutputVariables(variables{:});
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = WVModel(wvt)
            % Initialize a model from a WVTransform instance
            %
            % - Topic: Initialization
            % - Declaration: WVModel(wvt,options)
            % - Parameter wvt: a WaveVortexTranform instance
            %
            % 
            arguments
                wvt WVTransform {mustBeNonempty}
            end

            self.wvt = wvt; 
            self.initialOutputTime = self.t;
            self.initialTime = self.t;
            if self.wvt.hasClosure == false
                warning('The nonlinear flux has no damping and may not be stable.');
            end
            defaultGroup = self.addOutputGroup(self.defaultOutputGroupName);
            defaultGroup.setNetCDFOutputVariables(intersect({'Ap','Am','A0'},wvt.variableNames));
            self.particleIndexWithName = containers.Map();
            self.tracerIndexWithName = containers.Map();
        end
        
        function value = get.linearDynamics(self)
            value = isempty(self.wvt.forcing);
        end

        function value = get.t(self)
            value = self.wvt.t;
        end
        
        function addParticles(self,name,fluxOp,x,y,z,trackedFieldNames,options)
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
            arguments
                self WVModel {mustBeNonempty}
                name char {mustBeNonempty}
                fluxOp WVParticleFluxOperation {mustBeNonempty}
                x (1,:) double
                y (1,:) double
                z (1,:) double
            end
            arguments (Repeating)
                trackedFieldNames char
            end
            arguments
                options.trackedVarInterpolation char {mustBeMember(options.trackedVarInterpolation,["linear","spline","exact","finufft"])} = "spline"
                options.outputGroupName = "wave-vortex"
            end

            if self.didSetupIntegrator == 1
                error('Integrator was already setup, you cannot add particles.')
            end

            % Confirm that we really can track these variables.
            for iVar=1:length(trackedFieldNames)
                if ~any(ismember(self.wvt.variableNames,trackedFieldNames{iVar}))
                    error('Unable to find a WVVariableAnnotation named %s.', trackedFieldNames{iVar});
                end
                transformVar = self.wvt.propertyAnnotationWithName(trackedFieldNames{iVar});
                if isequal(self.wvt.spatialDimensionNames,{'x','y'})
                    if ~all(ismember(transformVar.dimensions,{'x','y'})) && ~all(ismember(transformVar.dimensions,{'x','y','z'}))
                        error('The WVVariableAnnotation %s does not have dimensions (x,y) or (x,y,z) and theforefore cannot be used for particle tracking', trackedFieldNames{iVar});
                    end
                else
                    if ~all(ismember(transformVar.dimensions,{'x','y','z'}))
                        error('The WVVariableAnnotation %s does not have dimensions x,y,z and theforefore cannot be used for particle tracking', trackedFieldNames{iVar});
                    end
                end
            end

            n = length(self.particle) + 1;

            self.particleIndexWithName(name) = n;
            self.particle{n}.name = name;
            self.particle{n}.fluxOp = fluxOp;
            self.particle{n}.x = x;
            self.particle{n}.y = y;
            self.particle{n}.z = z;
            self.particle{n}.trackedFieldNames = trackedFieldNames;
            trackedFields = struct;
            for i=1:length(trackedFieldNames)
                trackedFields.(trackedFieldNames{i}) = zeros(1,length(x));
            end
            self.particle{n}.trackedFields = trackedFields;
            self.particle{n}.trackedFieldInterpMethod = options.trackedVarInterpolation;

            self.updateParticleTrackedFields();

            outputGroup = self.outputGroupWithName(options.outputGroupName);
            outputGroup.addNetCDFOutputParticles(name);
        end

        function [x,y,z,trackedFields] = particlePositions(self,name)
            % Positions and values of tracked fields of particles at the current model time.
            %
            % - Topic: Particles
            % - Declaration: [x,y,z,trackedFields] = particlePositions(name)
            % - Parameter name: name of the particles
            p = self.particle{self.particleIndexWithName(name)};
            x = p.x;
            y = p.y;
            z = p.z;
            trackedFields = self.particle{self.particleIndexWithName(name)}.trackedFields;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Floats and drifters and tracer!
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
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
                options.outputGroupName = "wave-vortex"
            end
            floatFlux = WVParticleFluxOperation('floatFlux',@(wvt,x,y,z) wvt.variableAtPositionWithName(x,y,z,'u','v','w',interpolationMethod=options.advectionInterpolation));
            self.addParticles('float',floatFlux,x,y,z,trackedFields{:},trackedVarInterpolation=options.trackedVarInterpolation,outputGroupName=options.outputGroupName);
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
                options.outputGroupName = "wave-vortex"
            end
            drifterFlux = WVParticleFluxOperation('floatFlux',@(wvt,x,y,z) wvt.variableAtPositionWithName(x,y,z,'u','v',interpolationMethod=options.advectionInterpolation),isXYOnly=1);
            self.addParticles('drifter',drifterFlux,x,y,z,trackedFields{:},trackedVarInterpolation=options.trackedVarInterpolation,outputGroupName=options.outputGroupName);
        end

        function [x,y,z,tracked] = drifterPositions(self)
            % Current positions of the drifter particles
            % - Topic: Particles
            [x,y,z,tracked] = self.particlePositions('drifter');
        end

        function addTracer(self,phi,name)
            % Add a scalar field tracer to be advected by the flow
            % - Topic: Tracer
            n = length(self.tracerIndexWithName) + 1;
            self.tracerIndexWithName(name) = n;
            self.tracerArray{n} = phi;
        end

        function phi = tracer(self,name)
            % Scalar field of the requested tracer at the current model time.
            % - Topic: Tracer
            phi = self.tracerArray{self.tracerIndexWithName(name)};
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
            % - Parameter shouldUseScaledTolerance: (adapative) whether to scale by the energy norm. 0 or 1 (default)
            % - Parameter absToleranceA0: (adapative) absolute tolerance for A0 used when shouldUseScaledTolerance=0 . 1e-10 (default)
            % - Parameter absToleranceApm: (adapative) absolute tolerance for Apm used when shouldUseScaledTolerance=0 . 1e-10 (default)
            % - Parameter absToleranceXY: (adapative) absolute tolerance in meters for particle advection in (x,y). 1e-1 (default)
            % - Parameter absToleranceZ: (adapative) absolute tolerance  in meters for particle advection in (z). 1e-2 (default)
            % - Parameter shouldShowIntegrationStats: (adapative) whether to show integration output 0 or 1 (default)
            arguments
                self WVModel {mustBeNonempty}

                options.integratorType char {mustBeMember(options.integratorType,["fixed","adaptive"])} = "adaptive"

                fixedTimeStepOptions.deltaT (1,1) double {mustBePositive}
                fixedTimeStepOptions.cfl (1,1) double
                fixedTimeStepOptions.timeStepConstraint char {mustBeMember(fixedTimeStepOptions.timeStepConstraint,["advective","oscillatory","min"])} = "advective"

                adaptiveTimeStepOptions.integrator = @ode78
                adaptiveTimeStepOptions.absTolerance = 1e-6
                adaptiveTimeStepOptions.relTolerance = 1e-3;
                adaptiveTimeStepOptions.shouldUseScaledTolerance = 1;
                adaptiveTimeStepOptions.absToleranceA0 = 1e-10
                adaptiveTimeStepOptions.absToleranceApm = 1e-6
                adaptiveTimeStepOptions.absToleranceXY = 1e-1; % 100 km * 10^{-6}
                adaptiveTimeStepOptions.absToleranceZ = 1e-2;  
                adaptiveTimeStepOptions.shouldShowIntegrationStats double {mustBeMember(adaptiveTimeStepOptions.shouldShowIntegrationStats,[0 1])} = 0
            end

            self.resetFixedTimeStepIntegrator();
            self.resetAdapativeTimeStepIntegrator();

            self.integratorType = options.integratorType;
            if strcmp(self.integratorType,"adaptive")
                optionArgs = namedargs2cell(adaptiveTimeStepOptions);
                self.setupAdaptiveTimeStepIntegrator(optionArgs{:});
            else
                optionArgs = namedargs2cell(fixedTimeStepOptions);
                self.setupFixedTimeStepIntegrator(optionArgs{:});
            end 
        end
        
        function integrateToTime(self,finalTime,options)
            % Time step the model forward to the requested time.
            % - Topic: Integration
            arguments
                self WVModel {mustBeNonempty}
                finalTime (1,:) double
                options.shouldShowIntegrationDiagnostics double {mustBeMember(options.shouldShowIntegrationDiagnostics,[0 1])} = 1
            end
            if finalTime <= self.t
                fprintf('Reqested integration to time %d, but the model is currently at time t=%d.\n',round(finalTime),round(self.t));
                return;
            end
            if self.didSetupIntegrator ~= 1
                self.setupIntegrator();
            end
            self.shouldShowIntegrationDiagnostics = options.shouldShowIntegrationDiagnostics;
                  
            self.openNetCDFFileForTimeStepping();

            if strcmp(self.integratorType,"adaptive")
                self.integrateToTimeWithAdaptiveTimeStep(finalTime)
            else
                self.integrateToTimeWithFixedTimeStep(finalTime);
            end

            self.recordNetCDFFileHistory();            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % NetCDF Output
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function ncfile = createNetCDFFileForModelOutput(self,netcdfFile,options)
            % Create a NetCDF file for model output
            % - Topic: Writing to NetCDF files
            arguments
                self WVModel {mustBeNonempty}
                netcdfFile char {mustBeNonempty}
                options.outputInterval (1,1) double {mustBePositive}
                options.shouldOverwriteExisting logical = false
            end

            defaultGroup = self.outputGroupWithName(self.defaultOutputGroupName);
            if isempty(defaultGroup.outputInterval)
                if ~isfield(options,"outputInterval")
                    error("You must set an output interval");
                end
                defaultGroup.outputInterval = options.outputInterval;
            end

            hasOutputInterval = true;
            for iGroup =1:length(self.outputGroups)
                hasOutputInterval = hasOutputInterval & ~isempty(self.outputGroups(iGroup).outputInterval);
            end
            if ~hasOutputInterval
                error("One or more groups is missing an output interval");
            end

            properties = setdiff(self.wvt.requiredProperties,{'Ap','Am','A0','t'});

            ncfile = self.wvt.writeToFile(netcdfFile,properties{:},shouldOverwriteExisting=options.shouldOverwriteExisting,shouldAddRequiredProperties=false);
            self.ncfile = ncfile;
            self.didInitializeNetCDFFile = 0;
        end

        function recordNetCDFFileHistory(self,options)
            arguments
                self WVModel {mustBeNonempty}
                options.didBlowUp {mustBeNumeric} = 0
            end
            if isempty(self.ncfile)
                return
            end

            if options.didBlowUp == 1
                a = sprintf('%s: wrote %d time points to file. Terminated due to model blow-up.',datetime('now'),self.incrementsWrittenToFile);
            else
                a = sprintf('%s: wrote %d time points to file',datetime('now'),self.incrementsWrittenToFile);
            end
            if isKey(self.ncfile.attributes,'history')
                history = reshape(self.ncfile.attributes('history'),1,[]);
                history =cat(2,squeeze(history),a);
            else
                history = a;
            end
            self.ncfile.addAttribute('history',history);
        end

    end


    properties %(Access = protected)
        particle = {}  % cell array containing particle structs
        particleIndexWithName % map from particle name to cell array index

        tracerIndexWithName
        tracerArray = {}

        didSetupIntegrator=0
        didInitializeNetCDFFile=0

        integrationStartTime
        integrationLastInformWallTime       % wall clock, to keep track of the expected integration time
        integrationLastInformModelTime
        integrationInformTime = 10

        initialOutputTime   % output time corresponding to outputIndex=1 (set on instance initialization)

        integratorType      % Array integrator
        finalIntegrationTime % set only during an integration
        shouldShowIntegrationDiagnostics = 1

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

        function showIntegrationStartDiagnostics(self,finalTime)
            if self.shouldShowIntegrationDiagnostics  == 0
                return;
            end
            self.integrationStartTime = datetime('now');
            fprintf('Starting numerical simulation on %s.\n', datetime(self.integrationStartTime,TimeZone='local',Format='d-MMM-y HH:mm:ss Z'));
            fprintf('\tStarting at model time t=%.2f inertial periods and integrating to t=%.2f inertial periods.\n',self.t/self.wvt.inertialPeriod,finalTime/self.wvt.inertialPeriod);
            self.integrationLastInformWallTime = datetime('now');
            self.integrationLastInformModelTime = self.wvt.t;
        end

        function showIntegrationTimeDiagnostics(self,finalTime)
            if self.shouldShowIntegrationDiagnostics  == 0
                return;
            end
            deltaWallTime = datetime('now')-self.integrationLastInformWallTime;
            if ( seconds(deltaWallTime) > self.integrationInformTime)
                wallTimePerModelTime = deltaWallTime / (self.wvt.t - self.integrationLastInformModelTime);
                wallTimeRemaining = wallTimePerModelTime*(finalTime - self.wvt.t);
                fprintf('\tmodel time t=%.2f inertial periods. Estimated time to reach %.2f inertial periods is %s (%s)\n', self.t/self.wvt.inertialPeriod, finalTime/self.wvt.inertialPeriod, wallTimeRemaining, datetime(datetime('now')+wallTimeRemaining,TimeZone='local',Format='d-MMM-y HH:mm:ss Z')) ;
                self.wvt.summarizeEnergyContent();

                self.integrationLastInformWallTime = datetime('now');
                self.integrationLastInformModelTime = self.wvt.t;
            end
        end

        function showIntegrationFinishDiagnostics(self)
            if self.shouldShowIntegrationDiagnostics  == 0
                return;
            end
            integrationTotalTime = datetime('now')-self.integrationStartTime;
            fprintf('Finished after time %s.\n', integrationTotalTime);
        end

        function updateParticleTrackedFields(self)
            % One special thing we have to do is log the particle
            % tracked fields
            for iParticle=1:length(self.particle)
                trackedFieldNames = self.particle{iParticle}.trackedFieldNames;
                if ~isempty(trackedFieldNames)
                    varLagrangianValues = cell(1,length(trackedFieldNames));
                    p = self.particle{iParticle};
                    [varLagrangianValues{:}] = self.wvt.variableAtPositionWithName(p.x,p.y,p.z,trackedFieldNames{:},interpolationMethod=self.particle{iParticle}.trackedFieldInterpMethod);
                    for i=1:length(trackedFieldNames)
                        self.particle{iParticle}.trackedFields.(trackedFieldNames{i}) = varLagrangianValues{i};
                    end
                end
            end
        end

        function openNetCDFFileForTimeStepping(self)
            arguments (Input)
                self WVModel {mustBeNonempty}
            end
            if ~isempty(self.ncfile) && self.didInitializeNetCDFFile == 0
                for iGroup =1:length(self.outputGroups)
                    self.outputGroups(iGroup).openNetCDFFileForTimeStepping(self.ncfile);
                end
                self.didInitializeNetCDFFile = 1;
                self.incrementsWrittenToFile = 0;
                self.writeTimeStepToNetCDFFile();
            end
        
        end

        function writeTimeStepToNetCDFFile(self)
            if ( ~isempty(self.ncfile) && self.t > self.timeOfLastIncrementWrittenToFile )
                self.updateParticleTrackedFields();
                for iGroup =1:length(self.outputGroups)
                    self.outputGroups(iGroup).writeTimeStepToNetCDFFile(self.t);
                end
                self.timeOfLastIncrementWrittenToFile = self.t;
            end
        end

        function closeNetCDFFile(self)
            if ~isempty(self.ncfile)
                for iGroup =1:length(self.outputGroups)
                    self.outputGroups(iGroup).closeNetCDFFile();
                end
                self.ncfile.close();
            end
        end
    end


end