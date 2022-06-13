classdef WaveVortexModel < handle
    %WaveVortexModel Tools for integrating (time-stepping)
    %the WaveVortexModel.
    %
    %   inttool = WaveVortexModel(wvm,t) creates a new
    %   integration tool for the model. If a time (t) is given, the model
    %   coefficients are assumed represent that model time. Otherwise t=0.
    %
    %   inttool = WaveVortexModel(existingModelOutput,restartIndex,shouldDoubleResolution)
    %   opens existing NetCDF output from the WaveVortexModel and uses that
    %   for a restart. restartIndex is optional, defaults to Inf (last time
    %   point). shouldDoubleResolution is optional, defaults to 0.
    %
    %   'shouldOverwriteExisting', default 0
    %   'shouldDoubleResolution', 0 or 1
    %   'restartIndex', index in existingModelOutput to use as restart.

    % TODO, May 5th, 2022
    % - add multiple tracers, store ids and names in struct?
    % - remove tracer variance at aliased wavenumbers?
    % - add scalar option for floats and drifters, e.g., save density or pv
    % - linear dynamics should only save the coefficients once (actually
    %   option should be to only write initial conditions)
    % - need method of doing fancy stuff during the integration loop
    % - want to write float and drifter paths to memory
    % - Maybe a list of variables (as enums) that we want to write
    % - Definitely want to output physical variables some time
    % - openNetCDFFileForTimeStepping should report expected file size

    % AAGH. Getting myself twisted in knots over the right API.
    % while ( tool.integrateToTime(finalTime) )
    %   tool.incrementForward()
    %   tool.writeToFile()
    %
    %
    % Conflicts: once you've chosen an outputInterval, etc., you can't be
    % allowed to reset it. Same with output file.
    % Deal is though, writing to file or writing to memory should have the
    % same loop.
    %
    % Usages:
    % 1) init, 2) integrate to some final value, 3) integrate to another
    % 1) init, 2) set output interval, 3) int to final value, with stops
    % along the way at the output intervals.
    % init, set output interval to netcdf file, 
    %
    % To set the deltaT you *need* an outputInterval.
    % - Do you need the netcdf output file first? I don't think so
    % So,
    % Initialization (allowed once):
    % Option 1: WaveVortexModel(wvm,t)
    % Option 2: WaveVortexModel(existingModelOutput)
    % 
    % Setup the integrator (allowed once):
    % Option 1: <nothing>
    % Option 2: setupIntegrator(deltaT)
    % Option 3: setupIntegrator(deltaT,outputInterval)
    % Option 4: setupIntegrator(deltaT,outputInterval,finalTime) -- alt: SetupIntegratorForFixedIntegrationTime
    % return estimated time steps?
    %
    % Integrate (called repeatedly):
    % Option 1: modelTime = integrateOneTimeStep()
    % Option 2: modelTime = integrateToNextOutputTime()
    % Option 3: modelTime = integrateToTime(futureTime)
    %
    % showIntegrationTimeDiagnostics( ???? )
    %
    %
    % Setup NetCDF
    % createNetCDFFileForModelOutput
    % AppendToExisting

    properties
        wvt             % WaveVortexTransform
        t=0             % current model time (in seconds)
        initialTime=0
    
        nonlinearFlux      % TransformOperation (for now anyway).

        outputInterval      % model output interval (seconds)
        initialOutputTime   % output time corresponding to outputIndex=1 (set on instance initialization)

        integrator      % Array integrator
        
        % These methods all assume a fixed time-step integrator
        startTime       % wall clock, to keep track of the expected integration time
        stepsTaken=0    % number of RK4 steps/increments that have been made
        nSteps=inf      % total number of expected RK4 increments to reach the model time requested by the user

        stepsPerOutput      % number of RK4 steps between each output
        firstOutputStep     % first RK4 step that should be output. 0 indicates the initial conditions should be output
        outputIndex=1       % output index of the current/most recent step. If stepsTaken=0, outputIndex=1 means the initial conditions get written at index 1
        
        
        % *if* outputting to NetCDF file, these will be populated
        ncfile      % WaveVortexModelNetCDFFile instance---empty indicates no file output
        

        incrementsWrittenToFile
    end

    properties (SetAccess = private)
        particle = {}  % cell array containing particle structs
        particleIndexWithName % map from particle name to cell array index

        tracerIndexWithName
        tracerArray = {}

        didSetupIntegrator=0
        variablesToWriteToFile = {}
        initialConditionOnlyVariables = {}
        timeSeriesVariables = {}

        netcdfVariableMapForParticleWithName % map to a map containing the particle variables, e.g. particlesWithName('float') returns a map containing keys ('x','y','z') at minimum
    end

    properties (Dependent)
        linearDynamics
    end

    methods

        function WriteVariablesToFile(self,variables)
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
            self.variablesToWriteToFile = variables;
        end



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = WaveVortexModel(wvt)
            arguments
                wvt WaveVortexTransform {mustBeNonempty}
            end

            self.wvt = wvt; 
            self.t = wvt.t;
            self.initialOutputTime = self.t;
            self.initialTime = self.t;
             

%             nlFlux = NonlinearBoussinesqWithReducedInteractionMasks(self.wvt);
            self.nonlinearFlux = SingleModeQGPVE(self.wvt);
%             self.wvt.addTransformOperation(nlFlux);

            self.particleIndexWithName = containers.Map();
            self.tracerIndexWithName = containers.Map();
            self.netcdfVariableMapForParticleWithName = containers.Map();
        end
        

        function value = get.linearDynamics(self)
            value = isempty(self.nonlinearFlux);
        end
        
        function addParticles(self,name,fluxOp,x,y,z,trackedFieldNames,options)
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

            % Confirm that we really can track these variables.
            for iVar=1:length(trackedFieldNames)
                if ~isKey(self.wvt.stateVariableWithName,trackedFieldNames{iVar})
                    error('Unable to find a StateVariable named %s.', trackedFieldNames{iVar});
                end
                transformVar = self.wvt.stateVariableWithName(trackedFieldNames{iVar});
                if self.wvt.isBarotropic == 1
                    if ~all(ismember(transformVar.dimensions,{'x','y'})) && ~all(ismember(transformVar.dimensions,{'x','y','z'}))
                        error('The StateVariable %s does not have dimensions x,y and theforefore cannot be used for particle tracking', trackedFieldNames{iVar});
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


        function updateParticleTrackedFields(self)
            % One special thing we have to do is log the particle
            % tracked fields
            for iParticle=1:length(self.particle)
                trackedFieldNames = self.particle{iParticle}.trackedFieldNames;
                if ~isempty(trackedFieldNames)
                    varLagrangianValues = cell(1,length(trackedFieldNames));
                    p = self.particle{iParticle};
                    [varLagrangianValues{:}] = self.wvt.variablesAtPosition(p.x,p.y,p.z,trackedFieldNames{:},InterpolationMethod=self.particle{iParticle}.trackedFieldInterpMethod);
                    self.particle{iParticle}.trackedFields = vertcat(varLagrangianValues{:});
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Floats and drifters and tracer!
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function setFloatPositions(self,x,y,z,trackedFields,options)
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



        function stirWithConstituents(self,constituents)

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Integration loop
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function integrateToTime(self,finalTime,options)
            arguments
                self WaveVortexModel {mustBeNonempty}
                finalTime (1,:) double
                options.cfl (1,:) double = 0.5
                options.timeStepConstraint char {mustBeMember(options.timeStepConstraint,["advective","oscillatory","min"])} = "advective"
            end

            if self.didSetupIntegrator ~= 1
                [deltaT,advectiveDT,oscillatoryDT] = self.wvt.timeStepForCFL(options.cfl,self.outputInterval);
                if strcmp(options.timeStepConstraint,"advective")
                    deltaT = advectiveDT;
                elseif strcmp(options.timeStepConstraint,"oscillatory")
                    deltaT = oscillatoryDT;
                elseif strcmp(options.timeStepConstraint,"min")
                    deltaT = min(oscillatoryDT,advectiveDT);
                end
                self.setupIntegrator(deltaT,self.outputInterval,finalTime);
            end
   
            self.openNetCDFFileForTimeStepping();  
            while(self.t < finalTime)              
                self.integrateToNextOutputTime();
                self.writeTimeStepToNetCDFFile();
            end
        end

        function varargout = setupIntegrator(self,deltaT,outputInterval,finalTime)
            varargout = cell(1,0);
            if self.didSetupIntegrator == 1
                warning('You cannot setup the same integrator more than once.')
                return;
            end

            % logic through some default settings
            if nargin < 2 || isempty(deltaT) || deltaT <= 0
                didSetDeltaT = 0;
            else
                didSetDeltaT = 1;
            end
            
            if (nargin < 3 || isempty(outputInterval) || deltaT <= 0) && ~isempty(self.outputInterval)
                outputInterval = self.outputInterval;
            end
            if nargin < 3 || isempty(outputInterval)
                didSetOutputInterval = 0;
            else
                didSetOutputInterval = 1;
            end

            if didSetDeltaT == 0 && didSetOutputInterval == 0
                deltaT = self.wvt.timeStepForCFL(0.5);
            elseif didSetDeltaT == 0 && didSetOutputInterval == 1
                deltaT = self.wvt.timeStepForCFL(0.5,outputInterval);
            end
            
            if ~isempty(self.nonlinearFlux)
                self.wvt.addTransformOperation(self.nonlinearFlux);
            end

            % Now set the initial conditions and point the integrator to
            % the correct flux function
            Y0 = self.initialConditionsArray();
            if isempty(Y0{1})
                error('Nothing to do! You must have set to linear dynamics, without floats, drifters or tracers.');
            end
            self.integrator = ArrayIntegrator(@(t,y0) self.fluxAtTime(t,y0),Y0,deltaT);
            self.integrator.currentTime = self.t;

            if didSetOutputInterval == 1
                self.outputInterval = outputInterval;
                self.stepsPerOutput = round(outputInterval/self.integrator.stepSize);
                self.firstOutputStep = round((self.initialOutputTime-self.t)/self.integrator.stepSize);
            end

            if ~(nargin < 4 || isempty(finalTime))
                 % total dT time steps to meet or exceed the requested time.
                self.nSteps = ceil((finalTime-self.t)/self.integrator.stepSize);
                varargout{1} = length(self.firstOutputStep:self.stepsPerOutput:self.nSteps);
            end

            self.didSetupIntegrator = 1;
        end
        
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
            self.t = modelTime;
            if mod(self.stepsTaken - self.firstOutputStep,self.stepsPerOutput) == 0
                self.outputIndex = self.outputIndex + 1;

                self.updateParticleTrackedFields();
            end
        end


        function modelTime = integrateToNextOutputTime(self)
            if isempty(self.outputInterval)
                fprintf('You did not set an output interval, so how could I integrateToNextOutputTime?\n');
                return;
            end
            
            modelTime = self.integrateOneTimeStep;
            while( mod(self.stepsTaken - self.firstOutputStep,self.stepsPerOutput) ~= 0 )
                modelTime = self.integrateOneTimeStep;
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

        function showIntegrationTimeDiagnostics(self,integratorIncrement)
            if integratorIncrement == 0
                fprintf('Starting numerical simulation on %s.\n', datestr(datetime('now')));
                fprintf('\tStarting at model time t=%.2f inertial periods and integrating to t=%.2f inertial periods with %d RK4 time steps.\n',self.t/self.wvt.inertialPeriod,0/self.wvt.inertialPeriod,self.nSteps);
                if ~isempty(self.ncfile)
                    %fprintf('\tWriting %d of those time steps to file. Will write to output file starting at index %d.\n',sum(self.outputSteps>=0),self.outputIndex-self.incrementsWrittenToFile);
                end
            elseif integratorIncrement == 1
                self.startTime = datetime('now');
            else
                timePerStep = (datetime('now')-self.startTime)/(integratorIncrement-1);
                % We want to inform the user about every 30 seconds
                stepsPerInform = ceil(30/seconds(timePerStep));
                if (integratorIncrement==2 || mod(integratorIncrement,stepsPerInform) == 0)
                    timeRemaining = (self.nSteps-integratorIncrement+1)*timePerStep;
                    fprintf('\tmodel time t=%.2f inertial periods, RK4 time step %d of %d. Estimated finish time %s (%s from now)\n', self.t/inertialPeriod, integratorIncrement, self.nSteps, datestr(datetime('now')+timeRemaining), datestr(timeRemaining, 'HH:MM:SS')) ;
                    self.wvt.summarizeEnergyContent();
                end
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
            dx = (self.wvt.x(2)-self.wvt.x(1));

            advectiveDT = cfl*dx/U;
            oscillatoryDT = cfl*period;
            % A cfl of 1/12 for oscillatoryDT might be necessary for good numerical precision when advecting particles.

            fprintf('dX/U = %.1f s (%.1f min). The highest frequency resolved IGW has period of %.1f s (%.1f min).\n', dx/U,dx/U/60,period,period/60);

            if advectiveDT < oscillatoryDT
                deltaT = advectiveDT;
            else
                deltaT = oscillatoryDT;
            end

            if nargin == 3 && ~isempty(outputInterval)
                deltaT = outputInterval/ceil(outputInterval/deltaT);
                stepsPerOutput_ = round(outputInterval/deltaT);
                fprintf('Rounding to match the output interval dt: %.2f s (%d steps per output)\n',deltaT,stepsPerOutput_);
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

            ncfile = self.wvt.writeToFile(netcdfFile,shouldOverwriteExisting=options.shouldOverwriteExisting);

            % Now add a time dimension
            transformVar = self.wvt.stateVariableWithName('t');
            attributes = containers.Map();
            attributes('units') = transformVar.units;
            attributes('description') = transformVar.description;
            ncfile.addDimension(transformVar.name,[],attributes,options.Nt);

            if ~self.linearDynamics
                self.nonlinearFlux.writeToFile(ncfile,self.wvt);
            end

            self.ncfile = ncfile;
        end

        function openNetCDFFileForTimeStepping(self)
            arguments
                self WaveVortexModel {mustBeNonempty}
            end
%             if ~isempty(self.ncfile)
                % Sort through which variables we will record a time series
                % for, and which we will only write initial conditions.
                self.initialConditionOnlyVariables = {};
                self.timeSeriesVariables = {};
                for iVar = 1:length(self.variablesToWriteToFile)
                    transformVar = self.wvt.stateVariableWithName(self.variablesToWriteToFile{iVar});
                    attributes = containers.Map();
                    attributes('units') = transformVar.units;
                    attributes('description') = transformVar.description;

                    if (self.linearDynamics == 1 && transformVar.isVariableWithLinearTimeStep == 1) || (self.linearDynamics == 0 && transformVar.isVariableWithNonlinearTimeStep == 1)
                        self.timeSeriesVariables{end+1} = self.variablesToWriteToFile{iVar};
                        if transformVar.isComplex == 1
                            self.ncfile.initComplexVariable(transformVar.name,horzcat(transformVar.dimensions,'t'),attributes,'NC_DOUBLE');
                            self.ncfile.setVariable(transformVar.name,self.wvt.(transformVar.name));
                        else
                            self.ncfile.addVariable(transformVar.name,self.wvt.(transformVar.name),horzcat(transformVar.dimensions,'t'),attributes);
                        end
                    else
                        self.initialConditionOnlyVariables{end+1} = self.variablesToWriteToFile{iVar};
                        if transformVar.isComplex == 1
                            self.ncfile.initComplexVariable(transformVar.name,transformVar.dimensions,attributes,'NC_DOUBLE');
                            self.ncfile.setVariable(transformVar.name,self.wvt.(transformVar.name));
                        else
                            self.ncfile.addVariable(transformVar.name,self.wvt.(transformVar.name),transformVar.dimensions,attributes);
                        end
                    end   
                end

                for iTracer = 1:length(self.tracerArray)
                    if self.wvt.isBarotropic
                        self.ncfile.initVariable(self.tracerNames{iTracer}, {'x','y','t'},containers.Map({'isTracer'},{'1'}),'NC_DOUBLE');
                    else
                        self.ncfile.initVariable(self.tracerNames{iTracer}, {'x','y','z','t'},containers.Map({'isTracer'},{'1'}),'NC_DOUBLE');
                    end
                end

                for iParticle = 1:length(self.particle)
                    self.initializeParticleStorage(self.particle{iParticle}.name,size(self.particle{iParticle}.x,2),self.particle{iParticle}.trackedFieldNames{:});
                end

%             else
%                 if isempty(self.ncfile.ncid)
%                     self.ncfile.open();
%                 end
%             end

            self.incrementsWrittenToFile = 0;

            % Save the initial conditions
            self.writeTimeStepToNetCDFFile();         
        end

        function writeTimeStepToNetCDFFile(self)
            if ( ~isempty(self.ncfile) && mod(self.stepsTaken - self.firstOutputStep,self.stepsPerOutput) == 0 )
                self.ncfile.concatenateVariableAlongDimension('t',self.t,'t',self.outputIndex);

                for iVar=1:length(self.timeSeriesVariables)
                    selself.ncfilef.concatenateVariableAlongDimension(self.timeSeriesVariables{iVar},self.wvt.(self.timeSeriesVariables{iVar}),'t',self.outputIndex);
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
                    self.ncfile.concatenateVariableAlongDimension(strcat(particleName,'-',trackedFieldNames{iField}),trackedField.(trackedFieldNames{iField}),'t',iTime);
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Integration
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        

    end
end