classdef WaveVortexModel < handle
    %WaveVortexModel Tools for integrating (time-stepping)
    %the WaveVortexModel.
    %
    %   model = WaveVortexModel(wvt) creates a new model using the wvt
    %   (WaveVortexTransform)

    properties (GetAccess=public,SetAccess=protected)
        % The WaveVortexTransform instance the represents the ocean state.
        % Set on initialization only, the WaveVortexTransform in the model
        % performs all computations necessary to return information about
        % the ocean state at a given time.
        wvt WaveVortexTransform {mustBeNonempty}

        % Initial model time (seconds)
        % The time of the WaveVortexTransform when the model was
        % initialized. This also corresponds to the first time in the
        % NetCDF output file.
        initialTime (1,1) double = 0 
    
        % The operation responsible for computing the nonlinear flux of the model
        % If the nonlinearFlux is nil, then the model will advance using
        % linear dynamics (i.e., the wave-vortex coefficients will not
        % change).
        nonlinearFlux NonlinearFluxOperation
        
        % Reference to the NetCDFFile being used for model output
        % Empty indicates no file output.
        ncfile NetCDFFile
        
        % List of all StateVariables being written to NetCDF file
        % The default list includes 'Ap', 'Am', and 'A0' which is the
        % minimum set of variables required for a restart.
        netCDFOutputVariables = {'Ap','Am','A0'}

        % Model output interval (seconds)
        % This property is optionally set when calling setupIntegrator. If
        % set, it will allow you to call -integrateToNextOutputTime and, if
        % a NetCDF file is set for output, it will set the interval at
        % which time steps are written to file.
        outputInterval (1,1) double
    end

    properties (Dependent)
        % Indicates whether or not the model is using linear or nonlinear dynamics.
        % In practice, this is simply checking whether the nonlinearFlux
        % property is nil.
        linearDynamics

        % Current model time (seconds)
        % Current time of the ocean state, particle positions, and tracer.
        t (1,1) double
    end





    methods

        function addNetCDFOutputVariables(self,variables)
            % Add variables to list of variables to be written to the NetCDF variable during the model run.
            %
            % Pass strings of WaveVortexTransform state variables of the
            % same name. This must be called before using any of the
            % integrate methods.
            %
            % ```matlab
            % model.addNetCDFOutputVariables('A0','u','v');
            % ```
            %
            % - Parameter variables: strings of variable names.
            arguments
                self WaveVortexModel
            end
            arguments (Repeating)
                variables char
            end
            unknownVars = setdiff(variables,self.wvt.stateVariableWithName.keys);
            if ~isempty(unknownVars)
               error('The WaveVortexTransform does not have a variable named %s',unknownVars{1}) ;
            end
            self.netCDFOutputVariables = union(self.netCDFOutputVariables,variables);
        end

        function setNetCDFOutputVariables(self,variables)
            % Set list of variables to be written to the NetCDF variable during the model run.
            %
            % Pass strings of WaveVortexTransform state variables of the
            % same name. This must be called before using any of the
            % integrate methods.
            %
            % ```matlab
            % model.setNetCDFOutputVariables('A0','u','v');
            % ```
            %
            % - Parameter variables: strings of variable names.
            arguments
                self WaveVortexModel
            end
            arguments (Repeating)
                variables char
            end
            unknownVars = setdiff(variables,self.wvt.stateVariableWithName.keys);
            if ~isempty(unknownVars)
                error('The WaveVortexTransform does not have a variable named %s',unknownVars{1}) ;
            end
            self.netCDFOutputVariables = variables;
        end

        function removeNetCDFOutputVariables(self,variables)
            % Remove variables from the list of variables to be written to the NetCDF variable during the model run.
            %
            % Pass strings of WaveVortexTransform state variables of the
            % same name. This must be called before using any of the
            % integrate methods.
            %
            % ```matlab
            % model.removeNetCDFOutputVariables('A0','u','v');
            % ```
            %
            % - Parameter variables: strings of variable names.
            arguments
                self WaveVortexModel
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

        function self = WaveVortexModel(wvt,options)
            arguments
                wvt WaveVortexTransform {mustBeNonempty}
                options.nonlinearFlux NonlinearFluxOperation
            end

            self.wvt = wvt; 
            self.initialOutputTime = self.t;
            self.initialTime = self.t;
            if isfield(options,"nonlinearFlux")
                self.nonlinearFlux = options.nonlinearFlux;
                self.wvt.addTransformOperation(self.nonlinearFlux);
            end

%             nlFlux = NonlinearBoussinesqWithReducedInteractionMasks(self.wvt);
%             self.nonlinearFlux = SingleModeQGPVE(self.wvt);
%             self.wvt.addTransformOperation(nlFlux);

            self.particleIndexWithName = containers.Map();
            self.tracerIndexWithName = containers.Map();
            self.netcdfVariableMapForParticleWithName = containers.Map();
        end
        

        function value = get.linearDynamics(self)
            value = isempty(self.nonlinearFlux);
        end

        function value = get.t(self)
            value = self.wvt.t;
        end
        
        function addParticles(self,name,fluxOp,x,y,z,trackedFieldNames,options)
            % Add particles to be advected by the flow.
            arguments
                self WaveVortexModel {mustBeNonempty}
                name char {mustBeNonempty}
                fluxOp ParticleFluxOperation {mustBeNonempty}
                x (1,:) double
                y (1,:) double
                z (1,:) double
            end
            arguments (Repeating)
                trackedFieldNames char
            end
            arguments
                options.trackedVarInterpolation char {mustBeMember(options.trackedVarInterpolation,["linear","spline","exact"])} = "spline"
            end

            if self.didSetupIntegrator == 1
                error('Integrator was already setup, you cannot add particles.')
            end

            % Confirm that we really can track these variables.
            for iVar=1:length(trackedFieldNames)
                if ~isKey(self.wvt.stateVariableWithName,trackedFieldNames{iVar})
                    error('Unable to find a StateVariable named %s.', trackedFieldNames{iVar});
                end
                transformVar = self.wvt.stateVariableWithName(trackedFieldNames{iVar});
                if self.wvt.isBarotropic == 1
                    if ~all(ismember(transformVar.dimensions,{'x','y'})) && ~all(ismember(transformVar.dimensions,{'x','y','z'}))
                        error('The StateVariable %s does not have dimensions (x,y) or (x,y,z) and theforefore cannot be used for particle tracking', trackedFieldNames{iVar});
                    end
                else
                    if ~all(ismember(transformVar.dimensions,{'x','y','z'}))
                        error('The StateVariable %s does not have dimensions x,y,z and theforefore cannot be used for particle tracking', trackedFieldNames{iVar});
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
            % Set positions of float-like particles to be advected by the
            % model.
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
            % You can track the value of any known StateVariable along the
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
            %
            % - Parameter x: x-coordinate location of the particles
            % - Parameter y: y-coordinate location of the particles
            % - Parameter z: z-coordinate location of the particles
            % - Parameter trackedFields: strings of variable names
            % - Parameter advectionInterpolation: (optional) interpolation method used for particle advection. "linear" (default),"spline","exact"
            % - Parameter trackedVarInterpolation: (optional) interpolation method used for tracked field. "linear" (default),"spline","exact"
            arguments
                self WaveVortexModel {mustBeNonempty}
                x (1,:) double
                y (1,:) double
                z (1,:) double = []
            end
            arguments (Repeating)
                trackedFields char
            end
            arguments
                options.advectionInterpolation char {mustBeMember(options.advectionInterpolation,["linear","spline","exact"])} = "linear"
                options.trackedVarInterpolation char {mustBeMember(options.trackedVarInterpolation,["linear","spline","exact"])} = "linear"
            end
            floatFlux = ParticleFluxOperation('floatFlux',@(wvt,x,y,z) wvt.variablesAtPosition(x,y,z,'u','v','w',InterpolationMethod=options.advectionInterpolation));
            self.addParticles('float',floatFlux,x,y,z,trackedFields{:},trackedVarInterpolation=options.trackedVarInterpolation);
        end

        function [x,y,z,tracked] = floatPositions(self)
            % Returns the positions of the floats at the current time as well as the value of the fields being tracked.
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
            arguments
                self WaveVortexModel {mustBeNonempty}
                x (1,:) double
                y (1,:) double
                z (1,:) double = []
            end
            arguments (Repeating)
                trackedFields char
            end
            arguments
                options.advectionInterpolation char {mustBeMember(options.advectionInterpolation,["linear","spline","exact"])} = "linear"
                options.trackedVarInterpolation char {mustBeMember(options.trackedVarInterpolation,["linear","spline","exact"])} = "linear"
            end
            drifterFlux = ParticleFluxOperation('floatFlux',@(wvt,x,y,z) wvt.variablesAtPosition(x,y,z,'u','v',InterpolationMethod=options.advectionInterpolation),isXYOnly=1);
            self.addParticles('drifter',drifterFlux,x,y,z,trackedFields{:},trackedVarInterpolation=options.trackedVarInterpolation);
        end

        function [x,y,z,tracked] = drifterPositions(self)
            [x,y,z,tracked] = self.particlePositions('drifter');
        end

        function addTracer(self,phi,name)
            n = length(self.tracerIndexWithName) + 1;
            self.tracerIndexWithName(name) = n;
            self.tracerArray{n} = phi;
        end

        function phi = tracer(self,name)
            phi = self.tracerArray{self.tracerIndexWithName(name)};
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Integration loop
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function varargout = setupIntegrator(self,options)
            arguments
                self WaveVortexModel {mustBeNonempty}
                options.deltaT (1,1) double {mustBePositive}
                options.cfl (1,1) double
                options.timeStepConstraint char {mustBeMember(options.timeStepConstraint,["advective","oscillatory","min"])} = "advective"
                options.outputInterval
                options.finalTime
            end
            
            if self.didSetupIntegrator == 1
                warning('You cannot setup the same integrator more than once.')
                return;
            end

            if ~isempty(self.ncfile) && ~isfield(options,"outputInterval")
                error('You set a NetCDF output file, but have not set an outputInterval. Use -setupIntegrator with an outputInterval.');
            end

            if isfield(options,"deltaT")
                if isfield(options,"cfl")
                    warning('deltaT was already set, ignoring cfl')
                end
                deltaT = options.deltaT;

                if isfield(options,"outputInterval")
                    deltaT = options.outputInterval/ceil(options.outputInterval/deltaT);
                    if options.deltaT ~= deltaT
                        warning('deltaT changed from %f to %f to match the outputInterval.',options.deltaT,deltaT);
                    end
                end
            else
                if ~isfield(options,"cfl")
                    options.cfl = 0.25;
                end
                if isfield(options,"outputInterval")
                    [deltaT,advectiveDT,oscillatoryDT] = self.timeStepForCFL(options.cfl,options.outputInterval);
                else
                    [deltaT,advectiveDT,oscillatoryDT] = self.timeStepForCFL(options.cfl);
                end
                if strcmp(options.timeStepConstraint,"advective")
                    deltaT = advectiveDT;
                    fprintf('Using the advective dt')
                
                elseif strcmp(options.timeStepConstraint,"oscillatory")
                    deltaT = oscillatoryDT;
                    fprintf('Using the oscillatory dt')
                elseif strcmp(options.timeStepConstraint,"min")
                    deltaT = min(oscillatoryDT,advectiveDT);
                    fprintf('Using the min dt')
                end
                if isfield(options,"outputInterval")
                    fprintf(': %.2f s (%d steps per output)\n',deltaT,round(options.outputInterval/deltaT));
                else
                    fprintf(': %.2f s\n',deltaT);
                end
            end

            % Now set the initial conditions and point the integrator to
            % the correct flux function
            Y0 = self.initialConditionsArray();
            if isempty(Y0{1})
                error('Nothing to do! You must have set to linear dynamics, without floats, drifters or tracers.');
            end
            self.integrator = ArrayIntegrator(@(t,y0) self.fluxAtTime(t,y0),Y0,deltaT,currentTime=self.t);

            if isfield(options,"outputInterval")
                self.outputInterval = options.outputInterval;
                self.stepsPerOutput = round(self.outputInterval/self.integrator.stepSize);
                self.firstOutputStep = round((self.initialOutputTime-self.t)/self.integrator.stepSize);
            end

            if isfield(options,"finalTime")
                 % total dT time steps to meet or exceed the requested time.
                self.nSteps = ceil((options.finalTime-self.t)/self.integrator.stepSize);
                varargout{1} = length(self.firstOutputStep:self.stepsPerOutput:self.nSteps);
            else
                varargout = cell(1,0);
            end

            self.didSetupIntegrator = 1;
        end
        

        function integrateToTime(self,finalTime)
            arguments
                self WaveVortexModel {mustBeNonempty}
                finalTime (1,:) double
            end

            if self.didSetupIntegrator ~= 1
                self.setupIntegrator();
            end
            
            self.showIntegrationStartDiagnostics(finalTime);
            self.openNetCDFFileForTimeStepping();
            while(self.t < finalTime)
                self.integrateOneTimeStep();
                if self.didBlowUp == 1
                    return;
                end
                self.showIntegrationTimeDiagnostics(finalTime);
            end
        end


        function modelTime = integrateToNextOutputTime(self)
            if isempty(self.outputInterval)
                fprintf('You did not set an output interval, so how could I integrateToNextOutputTime?\n');
                return;
            end
            
            self.openNetCDFFileForTimeStepping();
            if mod(self.stepsTaken - self.firstOutputStep,self.stepsPerOutput) == 0
                modelTime = self.integrateOneTimeStep;
                if self.didBlowUp == 1
                    return;
                end
            end

            while( mod(self.stepsTaken - self.firstOutputStep,self.stepsPerOutput) ~= 0 )
                modelTime = self.integrateOneTimeStep;
                if self.didBlowUp == 1
                    return;
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % NetCDF Output
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function ncfile = createNetCDFFileForModelOutput(self,netcdfFile,options)
            arguments
                self WaveVortexModel {mustBeNonempty}
                netcdfFile char {mustBeNonempty}
                options.Nt (1,1) double {mustBePositive} = Inf
                options.shouldOverwriteExisting (1,1) {mustBeNumeric} = 0
            end

            ncfile = self.wvt.writeToFile(netcdfFile,shouldOverwriteExisting=options.shouldOverwriteExisting,shouldAddDefaultVariables=0);

            % Now add a time dimension
            transformVar = self.wvt.stateVariableWithName('t');
            attributes = containers.Map();
            attributes('units') = transformVar.units;
            attributes('description') = transformVar.description;
            ncfile.addDimension(transformVar.name,[],attributes,options.Nt);

            if ~self.linearDynamics
                ncfile.addAttribute('NonlinearFluxOperation',class(self.nonlinearFlux));
                self.nonlinearFlux.writeToFile(ncfile,self.wvt);
            end

            self.ncfile = ncfile;
        end

    end


    properties (Access = protected)
        particle = {}  % cell array containing particle structs
        particleIndexWithName % map from particle name to cell array index

        tracerIndexWithName
        tracerArray = {}

        didSetupIntegrator=0
        didInitializeNetCDFFile=0
        
        initialConditionOnlyVariables = {}
        timeSeriesVariables = {}

        netcdfVariableMapForParticleWithName % map to a map containing the particle variables, e.g. particlesWithName('float') returns a map containing keys ('x','y','z') at minimum

        integrationLastInformWallTime       % wall clock, to keep track of the expected integration time
        integrationLastInformModelTime
        integrationInformTime = 10

        initialOutputTime   % output time corresponding to outputIndex=1 (set on instance initialization)

        integrator      % Array integrator

        % These methods all assume a fixed time-step integrator

        stepsTaken=0    % number of RK4 steps/increments that have been made
        nSteps=inf      % total number of expected RK4 increments to reach the model time requested by the user

        stepsPerOutput      % number of RK4 steps between each output
        firstOutputStep     % first RK4 step that should be output. 0 indicates the initial conditions should be output
        outputIndex=1       % output index of the current/most recent step. If stepsTaken=0, outputIndex=1 means the initial conditions get written at index 1

        incrementsWrittenToFile
    end
    

    methods (Access=protected)

        function flag = didBlowUp(self)
            if ( any(isnan(self.wvt.Ap)|isnan(self.wvt.Am)|isnan(self.wvt.A0)) )
                flag = 1;
                fprintf('Blowup detected. Aborting.');
            else
                flag = 0;
            end
        end

        function showIntegrationStartDiagnostics(self,finalTime)
            fprintf('Starting numerical simulation on %s.\n', datestr(datetime('now')));
            fprintf('\tStarting at model time t=%.2f inertial periods and integrating to t=%.2f inertial periods.\n',self.t/self.wvt.inertialPeriod,finalTime/self.wvt.inertialPeriod);
            self.integrationLastInformWallTime = datetime('now');
            self.integrationLastInformModelTime = self.wvt.t;
        end

        function showIntegrationTimeDiagnostics(self,finalTime)
            deltaWallTime = datetime('now')-self.integrationLastInformWallTime;
            if ( seconds(deltaWallTime) > self.integrationInformTime)
                wallTimePerModelTime = deltaWallTime / (self.wvt.t - self.integrationLastInformModelTime);
                wallTimeRemaining = wallTimePerModelTime*(finalTime - self.wvt.t);
                fprintf('\tmodel time t=%.2f inertial periods. Estimated time to reach %.2f inertial periods is %s (%s)\n', self.t/self.wvt.inertialPeriod, finalTime/self.wvt.inertialPeriod, datestr(wallTimeRemaining, 'HH:MM:SS'), datestr(datetime('now')+wallTimeRemaining)) ;
                self.wvt.summarizeEnergyContent();

A02 = (self.wvt.A0_TE_factor/self.wvt.h) .* (self.wvt.A0.*conj(self.wvt.A0));
u_rms_alt = sqrt(2*sum(A02(:)));

fprintf('***temp hack***: u_rms: %f\n',u_rms_alt);

                self.integrationLastInformWallTime = datetime('now');
                self.integrationLastInformModelTime = self.wvt.t;
            end
        end

        function [deltaT,advectiveDT,oscillatoryDT] = timeStepForCFL(self, cfl, outputInterval)
            % Return the time step (in seconds) to maintain the given cfl condition.
            % If the cfl condition is not given, 0.25 will be assumed.
            % If outputInterval is given, the time step will be rounded to evenly
            % divide the outputInterval.
            if nargin == 1
                cfl = 0.25;
            end

            omega = self.wvt.Omega;
            period = 2*pi/max(abs(omega(:)));
            [u,v] = self.wvt.velocityField();
            U = max(max(max( sqrt(u.*u + v.*v) )));
            dx = (3/2)*(self.wvt.x(2)-self.wvt.x(1));

            advectiveDT = cfl*dx/U;
            oscillatoryDT = cfl*period;
            % A cfl of 1/12 for oscillatoryDT might be necessary for good numerical precision when advecting particles.

            fprintf('dX/U = %.1f s (%.1f min). The highest frequency resolved IGW has period of %.1f s (%.1f min).\n', dx/U,dx/U/60,period,period/60);

            if nargin == 3 && ~isempty(outputInterval)
                advectiveDT = outputInterval/ceil(outputInterval/advectiveDT);
                oscillatoryDT = outputInterval/ceil(outputInterval/oscillatoryDT);
            end

            if advectiveDT < oscillatoryDT
                deltaT = advectiveDT;
            else
                deltaT = oscillatoryDT;
            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Integration: initial conditions, flux, and one time step
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function Y0 = initialConditionsArray(self)
            Y0 = cell(1,1);
            n = 0;
            if self.linearDynamics == 0
                if self.nonlinearFlux.doesFluxAp == 1
                    n=n+1;Y0{n} = self.wvt.Ap;
                end
                if self.nonlinearFlux.doesFluxAm == 1
                    n=n+1;Y0{n} = self.wvt.Am;
                end
                if self.nonlinearFlux.doesFluxA0 == 1
                    n=n+1;Y0{n} = self.wvt.A0;
                end
            end

            for iParticles=1:length(self.particle)
                p = self.particle{iParticles};
                n=n+1;Y0{n} = p.x;
                n=n+1;Y0{n} = p.y;
                if ~self.particle{iParticles}.fluxOp.isXYOnly
                    n=n+1;Y0{n} = p.z;
                end
            end

            for i=1:length(self.tracerArray)
                n=n+1;Y0{n} = self.tracerArray{i};
            end
        end

        function F = fluxAtTime(self,t,y0)
            F = cell(1,1);
            n = 0;
            self.wvt.t = t;
            if self.linearDynamics == 0
                nlF = cell(1,self.nonlinearFlux.nVarOut);
                [nlF{:}] = self.nonlinearFlux.Compute(self.wvt);
                if self.nonlinearFlux.doesFluxAp == 1
                    n=n+1; F{n} = nlF{n};
                end
                if self.nonlinearFlux.doesFluxAm == 1
                    n=n+1; F{n} = nlF{n};
                end
                if self.nonlinearFlux.doesFluxA0 == 1
                    n=n+1; F{n} = nlF{n};
                end
            else

            end

            for iParticles=1:length(self.particle)
                p = self.particle{iParticles};
                if self.particle{iParticles}.fluxOp.isXYOnly
                    [F{n+1},F{n+2}] = self.particle{iParticles}.fluxOp.Compute(self.wvt,p.x,p.y,p.z);
                    n=n+2;
                else
                    [F{n+1},F{n+2},F{n+3}] = self.particle{iParticles}.fluxOp.Compute(self.wvt,p.x,p.y,p.z);
                    n=n+3;
                end
            end

            if ~isempty(self.tracerArray)
                for i=1:length(self.tracerArray)
                    phibar = self.wvt.transformFromSpatialDomainWithF(y0{n+1});
                    [~,Phi_x,Phi_y,Phi_z] = self.wvt.transformToSpatialDomainWithFAllDerivatives(phibar);
                    n=n+1;F{n} = -self.wvt.u .* Phi_x - self.wvt.v.*Phi_y - self.wvt.w.*Phi_z;
                end
            end
        end

        function modelTime = integrateOneTimeStep(self)
            % Ask the integrator to take one step forward, then record the
            % results.
            self.integrator.IncrementForward();
            n=0;
            if self.linearDynamics == 0
                if self.nonlinearFlux.doesFluxAp == 1
                    n=n+1; self.wvt.Ap = self.integrator.currentY{n};
                end
                if self.nonlinearFlux.doesFluxAm == 1
                    n=n+1; self.wvt.Am = self.integrator.currentY{n};
                end
                if self.nonlinearFlux.doesFluxA0 == 1
                    n=n+1; self.wvt.A0 = self.integrator.currentY{n};
                end
            end

            for iParticles=1:length(self.particle)
                n=n+1; self.particle{iParticles}.x = self.integrator.currentY{n};
                n=n+1; self.particle{iParticles}.y = self.integrator.currentY{n};
                if ~self.particle{iParticles}.fluxOp.isXYOnly
                    n=n+1; self.particle{iParticles}.z = self.integrator.currentY{n};
                end
            end

            for iTracer=1:length(self.tracerArray)
                n=n+1; self.tracerArray{iTracer} = self.integrator.currentY{n};
            end
            

            % Rather than use the integrator time, which add floating point
            % numbers each time step, we multiple the steps taken by the
            % step size. This reduces rounding errors.
            self.stepsTaken = self.stepsTaken + 1;
            modelTime = self.initialTime + self.stepsTaken * self.integrator.stepSize;
            self.wvt.t = modelTime;
            if mod(self.stepsTaken - self.firstOutputStep,self.stepsPerOutput) == 0
                self.outputIndex = self.outputIndex + 1;

                self.updateParticleTrackedFields();
            end

            self.writeTimeStepToNetCDFFile(); 
        end


        function updateParticleTrackedFields(self)
            % One special thing we have to do is log the particle
            % tracked fields
            for iParticle=1:length(self.particle)
                trackedFieldNames = self.particle{iParticle}.trackedFieldNames;
                if ~isempty(trackedFieldNames)
                    varLagrangianValues = cell(1,length(trackedFieldNames));
                    p = self.particle{iParticle};
                    [varLagrangianValues{:}] = self.wvt.variablesAtPosition(p.x,p.y,p.z,trackedFieldNames{:},InterpolationMethod=self.particle{iParticle}.trackedFieldInterpMethod);
                    for i=1:length(trackedFieldNames)
                        self.particle{iParticle}.trackedFields.(trackedFieldNames{i}) = varLagrangianValues{i};
                    end
                end
            end
        end



        function openNetCDFFileForTimeStepping(self)
            arguments
                self WaveVortexModel {mustBeNonempty}
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

                    transformVar = self.wvt.stateVariableWithName(self.netCDFOutputVariables{iVar});
                    attributes = containers.Map();
                    attributes('units') = transformVar.units;
                    attributes('description') = transformVar.description;

                    if (self.linearDynamics == 1 && transformVar.isVariableWithLinearTimeStep == 1) || (self.linearDynamics == 0 && transformVar.isVariableWithNonlinearTimeStep == 1)
                        self.timeSeriesVariables{end+1} = self.netCDFOutputVariables{iVar};
                        if transformVar.isComplex == 1
                            self.ncfile.initComplexVariable(transformVar.name,horzcat(transformVar.dimensions,'t'),attributes,'NC_DOUBLE');
                        else
                            self.ncfile.initVariable(transformVar.name,horzcat(transformVar.dimensions,'t'),attributes,'NC_DOUBLE');
                        end
                    else
                        self.initialConditionOnlyVariables{end+1} = self.netCDFOutputVariables{iVar};
                        if transformVar.isComplex == 1
                            self.ncfile.initComplexVariable(transformVar.name,transformVar.dimensions,attributes,'NC_DOUBLE');
                            self.ncfile.setVariable(transformVar.name,self.wvt.(transformVar.name));
                        else
                            self.ncfile.addVariable(transformVar.name,self.wvt.(transformVar.name),transformVar.dimensions,attributes);
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

                % Save the initial conditions
                self.writeTimeStepToNetCDFFile();
            end
        
        end

        function writeTimeStepToNetCDFFile(self)
            if ( ~isempty(self.ncfile) && mod(self.stepsTaken - self.firstOutputStep,self.stepsPerOutput) == 0 )
                self.ncfile.concatenateVariableAlongDimension('t',self.t,'t',self.outputIndex);

                for iVar=1:length(self.timeSeriesVariables)
                    self.ncfile.concatenateVariableAlongDimension(self.timeSeriesVariables{iVar},self.wvt.(self.timeSeriesVariables{iVar}),'t',self.outputIndex);
                end

                for iParticle = 1:length(self.particle)
                    [x,y,z,trackedFields] = self.particlePositions(self.particle{iParticle}.name);
                    self.writeParticleDataAtTimeIndex(self.particle{iParticle}.name,self.outputIndex,x,y,z,trackedFields);
                end

                for iTracer = 1:length(self.tracerArray)
                    self.ncfile.WriteTracerWithNameTimeAtIndex(self.outputIndex,self.tracerNames{iTracer},self.tracerArray{iTracer});
                end

                self.incrementsWrittenToFile = self.incrementsWrittenToFile + 1;
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
                self WaveVortexModel
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
            [dim,var] = self.ncfile.addDimension(strcat(particleName,'-id'),(1:nParticles).',attributes);
            variables('id') = var;
            
            % careful to create a new object each time we init
            if self.wvt.isBarotropic == 1
                dimVars = {'x','y'};
            else
                dimVars = {'x','y','z'};
            end
            for iVar=1:length(dimVars)
                attributes = containers.Map(commonKeys,commonVals);
                attributes('units') = self.wvt.transformDimensionWithName(dimVars{iVar}).units;
                attributes('particleVariableName') = dimVars{iVar};
                variables(dimVars{iVar}) = self.ncfile.initVariable(strcat(particleName,'-',dimVars{iVar}),{dim.name,'t'},attributes,'NC_DOUBLE');
            end

            for iVar=1:length(trackedFieldNames)
                transformVar = self.wvt.stateVariableWithName(trackedFieldNames{iVar});
                attributes = containers.Map(commonKeys,commonVals);
                attributes('units') = transformVar.units;
                attributes('particleVariableName') = trackedFieldNames{iVar};
                variables(trackedFieldNames{iVar}) = self.ncfile.initVariable(strcat(particleName,'-',trackedFieldNames{iVar}),{dim.name,'t'},attributes,'NC_DOUBLE');
            end
 
            self.netcdfVariableMapForParticleWithName(particleName) = variables;
        end

        function writeParticleDataAtTimeIndex(self,particleName,iTime,x,y,z,trackedFields)
            self.ncfile.concatenateVariableAlongDimension(strcat(particleName,'-x'),x,'t',iTime);
            self.ncfile.concatenateVariableAlongDimension(strcat(particleName,'-y'),y,'t',iTime);
            if ~self.wvt.isBarotropic
                self.ncfile.concatenateVariableAlongDimension(strcat(particleName,'-z'),z,'t',iTime);
            end

            if ~isempty(trackedFields)
                trackedFieldNames = fieldnames(trackedFields);
                for iField=1:length(trackedFieldNames)
                    self.ncfile.concatenateVariableAlongDimension(strcat(particleName,'-',trackedFieldNames{iField}),trackedFields.(trackedFieldNames{iField}),'t',iTime);
                end
            end
        end
    end

    methods (Static)
        model = modelFromFile(path,options)
    end
end