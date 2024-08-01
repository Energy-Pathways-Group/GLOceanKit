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
    % model = WVModel(wvt,nonlinearFlux=QGPVE(wvt,u_damp=wvt.uvMax));
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
    
    properties (GetAccess=public, SetAccess=public)
            % The operation responsible for computing the nonlinear flux of the model
            % - Topic: Model Properties
            % If the nonlinearFlux is nil, then the model will advance using
            % linear dynamics (i.e., the wave-vortex coefficients will not
            % change).
            nonlinearFluxOperation WVNonlinearFluxOperation
    end

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
        
        % List of all StateVariables being written to NetCDF file
        % - Topic: Writing to NetCDF files
        % The default list includes 'Ap', 'Am', and 'A0' which is the
        % minimum set of variables required for a restart.
        netCDFOutputVariables = {'Ap','Am','A0'}

        % Model output interval (seconds)
        % - Topic: Integration
        % This property is optionally set when calling setupIntegrator. If
        % set, it will allow you to call -integrateToNextOutputTime and, if
        % a NetCDF file is set for output, it will set the interval at
        % which time steps are written to file.
        outputInterval = [] % (1,1) double

        % output index of the current/most recent step.
        % - Topic: Integration
        % If stepsTaken=0, outputIndex=1 means the initial conditions get written at index 1
        incrementsWrittenToFile (1,1) uint64 = 0

        % output index of the current/most recent step.
        % - Topic: Integration
        % If stepsTaken=0, outputIndex=1 means the initial conditions get written at index 1
        timeOfLastIncrementWrittenToFile (1,1) double = Inf

        % output index of the anticipated next step.
        % - Topic: Integration
        % If stepsTaken=0, outputIndex=1 means the initial conditions get written at index 1
        timeOfNextIncrementToWriteToFile (1,1) double = Inf
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
    end


    methods (Static)
        model = modelFromFile(path,options)
    end

    methods

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
            unknownVars = setdiff(variables,self.wvt.variableNames);
            if ~isempty(unknownVars)
               error('The WVTransform does not have a variable named %s',unknownVars{1}) ;
            end
            self.netCDFOutputVariables = union(self.netCDFOutputVariables,variables);
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
            unknownVars = setdiff(variables,self.wvt.variableNames);
            if ~isempty(unknownVars)
                error('The WVTransform does not have a variable named %s',unknownVars{1}) ;
            end
            self.netCDFOutputVariables = variables;
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
            self.netCDFOutputVariables = setdiff(self.netCDFOutputVariables,variables);
        end

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
            % - Parameter nonlinearFlux: (optional) a WVNonlinearFluxOperation used to time-step the WVTransform forward in time.
            %
            % 
            arguments
                wvt WVTransform {mustBeNonempty}
                options.nonlinearFlux WVNonlinearFluxOperation
            end

            self.wvt = wvt; 
            self.initialOutputTime = self.t;
            self.initialTime = self.t;
            if isfield(options,"nonlinearFlux")
                self.nonlinearFluxOperation = options.nonlinearFlux;
                self.wvt.nonlinearFluxOperation = options.nonlinearFlux;
            else
                self.nonlinearFluxOperation = self.wvt.nonlinearFluxOperation;
                if self.nonlinearFluxOperation.nu_xy == 0
                    warning('The nonlinear flux has no damping.');
                end
            end

            self.particleIndexWithName = containers.Map();
            self.tracerIndexWithName = containers.Map();
            self.netcdfVariableMapForParticleWithName = containers.Map();
        end
        

        function value = get.linearDynamics(self)
            value = isempty(self.nonlinearFluxOperation);
        end

        function set.nonlinearFluxOperation(self,value)
            if (self.didInitializeNetCDFFile == 1 || self.didSetupIntegrator == 1)
                error('You cannot change the nonlinearFlux after the integrator has been setup or the NetCDFFile has been created.');
            end
            self.nonlinearFluxOperation = value;
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
            end

            if self.didSetupIntegrator == 1
                error('Integrator was already setup, you cannot add particles.')
            end

            % Confirm that we really can track these variables.
            for iVar=1:length(trackedFieldNames)
                if ~any(ismember(self.wvt.variableNames,trackedFieldNames{iVar}))
                    error('Unable to find a WVVariableAnnotation named %s.', trackedFieldNames{iVar});
                end
                transformVar = self.wvt.variableAnnotationWithName(trackedFieldNames{iVar});
                if self.wvt.isBarotropic == 1
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
            end
            floatFlux = WVParticleFluxOperation('floatFlux',@(wvt,x,y,z) wvt.variablesAtPosition(x,y,z,'u','v','w',interpolationMethod=options.advectionInterpolation));
            self.addParticles('float',floatFlux,x,y,z,trackedFields{:},trackedVarInterpolation=options.trackedVarInterpolation);
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
            end
            drifterFlux = WVParticleFluxOperation('floatFlux',@(wvt,x,y,z) wvt.variablesAtPosition(x,y,z,'u','v',interpolationMethod=options.advectionInterpolation),isXYOnly=1);
            self.addParticles('drifter',drifterFlux,x,y,z,trackedFields{:},trackedVarInterpolation=options.trackedVarInterpolation);
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
            % - Topic: Integration
            % - Declaration: setupIntegrator(self,options)
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
            if ~isempty(self.ncfile)
                outputTimes = self.timeOfLastIncrementWrittenToFile:self.outputInterval:finalTime;
                if isscalar(outputTimes)
                   outputTimes = cat(2,outputTimes,finalTime);
                end
            else
                outputTimes = [self.t finalTime];
            end

            self.finalIntegrationTime = finalTime;
            if strcmp(self.integratorType,"adaptive")
                self.integrateToTimeWithAdaptiveTimeStep(outputTimes)
            else
                self.integrateToTimeWithFixedTimeStep(outputTimes);
            end
            self.finalIntegrationTime = [];

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
                options.outputInterval (1,1) double {mustBePositive} = 86400
                options.Nt (1,1) double {mustBePositive} = Inf
                options.shouldOverwriteExisting double {mustBeMember(options.shouldOverwriteExisting,[0 1])} = 0
                options.shouldUseClassicNetCDF double {mustBeMember(options.shouldUseClassicNetCDF,[0 1])} = 1 
            end

            self.outputInterval = options.outputInterval;
            self.timeOfNextIncrementToWriteToFile = 0;
            ncfile = self.wvt.writeToFile(netcdfFile,shouldOverwriteExisting=options.shouldOverwriteExisting,shouldAddDefaultVariables=0,shouldUseClassicNetCDF=options.shouldUseClassicNetCDF);

            % Now add a time dimension
            varAnnotation = self.wvt.variableAnnotationWithName('t');
            varAnnotation.attributes('units') = varAnnotation.units;
            varAnnotation.attributes('long_name') = varAnnotation.description;
            varAnnotation.attributes('standard_name') = 'time';
            varAnnotation.attributes('long_name') = 'time';
            varAnnotation.attributes('units') = 'seconds since 1970-01-01 00:00:00';
            varAnnotation.attributes('axis') = 'T';
            varAnnotation.attributes('calendar') = 'standard';
            ncfile.addDimension(varAnnotation.name,[],varAnnotation.attributes,options.Nt);

            ncfile.addAttribute('shouldUseLinearDynamics',uint8(self.linearDynamics));

            self.ncfile = ncfile;
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
                a = sprintf('%s: wrote %d time points to file. Terminated to do model blow-up.',datetime('now'),self.incrementsWrittenToFile);
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
        
        initialConditionOnlyVariables = {}
        timeSeriesVariables = {}

        netcdfVariableMapForParticleWithName % map to a map containing the particle variables, e.g. particlesWithName('float') returns a map containing keys ('x','y','z') at minimum

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

% A02 = (self.wvt.A0_TE_factor/self.wvt.h) .* (self.wvt.A0.*conj(self.wvt.A0));
% u_rms_alt = sqrt(2*sum(A02(:)));
% 
% fprintf('***temp hack***: u_rms: %f\n',u_rms_alt);

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
                    [varLagrangianValues{:}] = self.wvt.variablesAtPosition(p.x,p.y,p.z,trackedFieldNames{:},interpolationMethod=self.particle{iParticle}.trackedFieldInterpMethod);
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
                % Sort through which variables we will record a time series
                % for, and which we will only write initial conditions.
                self.initialConditionOnlyVariables = {};
                self.timeSeriesVariables = {};
                for iVar = 1:length(self.netCDFOutputVariables)
                    if isKey(self.ncfile.variableWithName,self.netCDFOutputVariables{iVar}) || isKey(self.ncfile.complexVariableWithName,self.netCDFOutputVariables{iVar})
                        continue;
                    end

                    varAnnotation = self.wvt.variableAnnotationWithName(self.netCDFOutputVariables{iVar});
                    varAnnotation.attributes('units') = varAnnotation.units;
                    varAnnotation.attributes('long_name') = varAnnotation.description;

                    if (self.linearDynamics == 1 && varAnnotation.isVariableWithLinearTimeStep == 1) || (self.linearDynamics == 0 && varAnnotation.isVariableWithNonlinearTimeStep == 1)
                        self.timeSeriesVariables{end+1} = self.netCDFOutputVariables{iVar};
                        if varAnnotation.isComplex == 1
                            self.ncfile.initComplexVariable(varAnnotation.name,horzcat(varAnnotation.dimensions,'t'),varAnnotation.attributes,'NC_DOUBLE');
                        else
                            self.ncfile.initVariable(varAnnotation.name,horzcat(varAnnotation.dimensions,'t'),varAnnotation.attributes,'NC_DOUBLE');
                        end
                    else
                        self.initialConditionOnlyVariables{end+1} = self.netCDFOutputVariables{iVar};
                        if varAnnotation.isComplex == 1
                            self.ncfile.initComplexVariable(varAnnotation.name,varAnnotation.dimensions,varAnnotation.attributes,'NC_DOUBLE');
                            self.ncfile.setVariable(varAnnotation.name,self.wvt.(varAnnotation.name));
                        else
                            self.ncfile.addVariable(varAnnotation.name,self.wvt.(varAnnotation.name),varAnnotation.dimensions,varAnnotation.attributes);
                        end
                    end
                end

                for iTracer = 1:length(self.tracerArray)
                    if isKey(self.ncfile.variableWithName,self.tracerNames{iTracer})
                        continue;
                    end
                    if self.wvt.isBarotropic
                        self.ncfile.initVariable(self.tracerNames{iTracer}, {'x','y','t'},containers.Map({'isTracer'},{'1'}),'NC_DOUBLE');
                    else
                        self.ncfile.initVariable(self.tracerNames{iTracer}, {'x','y','z','t'},containers.Map({'isTracer'},{'1'}),'NC_DOUBLE');
                    end
                end

                for iParticle = 1:length(self.particle)
                    if isKey(self.ncfile.variableWithName,self.particle{iParticle}.name)
                        continue;
                    end
                    self.initializeParticleStorage(self.particle{iParticle}.name,size(self.particle{iParticle}.x,2),self.particle{iParticle}.trackedFieldNames{:});
                end

                self.didInitializeNetCDFFile = 1;
                self.incrementsWrittenToFile = 0;

                self.writeTimeStepToNetCDFFile();
            end
        
        end

        function writeTimeStepToNetCDFFile(self)
            if ( ~isempty(self.ncfile) && self.timeOfNextIncrementToWriteToFile == self.t )
                self.updateParticleTrackedFields();

                outputIndex = self.incrementsWrittenToFile + 1;

                self.ncfile.concatenateVariableAlongDimension('t',self.t,'t',outputIndex);

                for iVar=1:length(self.timeSeriesVariables)
                    self.ncfile.concatenateVariableAlongDimension(self.timeSeriesVariables{iVar},self.wvt.(self.timeSeriesVariables{iVar}),'t',outputIndex);
                end

                for iParticle = 1:length(self.particle)
                    [x,y,z,trackedFields] = self.particlePositions(self.particle{iParticle}.name);
                    self.writeParticleDataAtTimeIndex(self.particle{iParticle}.name,outputIndex,x,y,z,trackedFields);
                end

                for iTracer = 1:length(self.tracerArray)
                    self.ncfile.WriteTracerWithNameTimeAtIndex(outputIndex,self.tracerNames{iTracer},self.tracerArray{iTracer});
                end

                self.incrementsWrittenToFile = outputIndex;
                self.timeOfLastIncrementWrittenToFile = self.t;
                self.timeOfNextIncrementToWriteToFile = self.t + self.outputInterval;
            end
        end

        function closeNetCDFFile(self)
            if ~isempty(self.ncfile)
                fprintf('Ending simulation. Wrote %d time points to file\n',self.incrementsWrittenToFile);
                self.ncfile.close();
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Particles
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function initializeParticleStorage(self,particleName, nParticles, trackedFieldNames)
            arguments
                self WVModel
                particleName char
                nParticles (1,1) double {mustBePositive}
            end
            arguments (Repeating)
                trackedFieldNames char
            end

            variables = containers.Map();

            commonKeys = {'isParticle','particleName'};
            commonVals = {1,particleName};
            attributes = containers.Map(commonKeys,commonVals);
            attributes('units') = 'unitless id number';
            attributes('particleVariableName') = 'id';
            [dim,var] = self.ncfile.addDimension(strcat(particleName,'_id'),(1:nParticles).',attributes);
            variables('id') = var;
            
            % careful to create a new object each time we init
            if self.wvt.isBarotropic == 1
                dimVars = {'x','y'};
            else
                dimVars = {'x','y','z'};
            end
            for iVar=1:length(dimVars)
                attributes = containers.Map(commonKeys,commonVals);
                attributes('units') = self.wvt.dimensionAnnotationWithName(dimVars{iVar}).units;
                attributes('long_name') = strcat(self.wvt.dimensionAnnotationWithName(dimVars{iVar}).description,', recorded along the particle trajectory');
                attributes('particleVariableName') = dimVars{iVar};
                variables(dimVars{iVar}) = self.ncfile.initVariable(strcat(particleName,'_',dimVars{iVar}),{dim.name,'t'},attributes,'NC_DOUBLE');
            end

            for iVar=1:length(trackedFieldNames)
                varAnnotation = self.wvt.variableAnnotationWithName(trackedFieldNames{iVar});
                attributes = containers.Map(commonKeys,commonVals);
                attributes('units') = varAnnotation.units;
                attributes('long_name') = strcat(varAnnotation.description,', recorded along the particle trajectory');
                attributes('particleVariableName') = trackedFieldNames{iVar};
                variables(trackedFieldNames{iVar}) = self.ncfile.initVariable(strcat(particleName,'_',trackedFieldNames{iVar}),{dim.name,'t'},attributes,'NC_DOUBLE');
            end
 
            self.netcdfVariableMapForParticleWithName(particleName) = variables;
        end

        function writeParticleDataAtTimeIndex(self,particleName,iTime,x,y,z,trackedFields)
            self.ncfile.concatenateVariableAlongDimension(strcat(particleName,'_x'),x,'t',iTime);
            self.ncfile.concatenateVariableAlongDimension(strcat(particleName,'_y'),y,'t',iTime);
            if ~self.wvt.isBarotropic
                self.ncfile.concatenateVariableAlongDimension(strcat(particleName,'_z'),z,'t',iTime);
            end

            if ~isempty(trackedFields)
                trackedFieldNames = fieldnames(trackedFields);
                for iField=1:length(trackedFieldNames)
                    self.ncfile.concatenateVariableAlongDimension(strcat(particleName,'_',trackedFieldNames{iField}),trackedFields.(trackedFieldNames{iField}),'t',iTime);
                end
            end
        end
    end


end