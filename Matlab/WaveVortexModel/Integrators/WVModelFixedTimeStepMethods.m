classdef WVModelFixedTimeStepMethods < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    % properties (Abstract,GetAccess=public, SetAccess=public)
    %     nonlinearFluxOperation
    % end
    properties (Abstract,GetAccess=public, SetAccess=protected)
        wvt
    end
    properties (Abstract) %(Access = protected)
        particle
        tracerArray
        linearDynamics
        t
        didSetupIntegrator
        finalIntegrationTime
        outputGroups
    end
    properties
        integratorOptions
        rk4Integrator
        % These methods all assume a fixed time-step integrator
        stepsTaken=0    % number of RK4 steps/increments that have been made
        nSteps=inf      % total number of expected RK4 increments to reach the model time requested by the user

        stepsPerOutput      % number of RK4 steps between each output
        firstOutputStep     % first RK4 step that should be output. 0 indicates the initial conditions should be output

        % Model output interval (seconds)
        % - Topic: Integration
        % This property is optionally set when calling setupIntegrator. If
        % set, it will allow you to call -integrateToNextOutputTime and, if
        % a NetCDF file is set for output, it will set the interval at
        % which time steps are written to file.
        outputInterval
    end

    methods (Abstract)
        flag = didBlowUp(self)
        showIntegrationStartDiagnostics(self,finalTime)
        showIntegrationTimeDiagnostics(self,finalTime)
        showIntegrationFinishDiagnostics(self)
        outputGroupWithName(self,name)
    end
    methods (Static,Abstract)
        name = defaultOutputGroupName()
    end

    methods

        function setupFixedTimeStepIntegrator(self,options)
            % Customize the time-stepping
            %
            % - Topic: Integration
            % - Declaration: setupIntegrator(self,options)
            % - Parameter deltaT: (optional) set the integrator time step
            % - Parameter cfl: (optional) set the cfl condition used to set integrator time step
            % - Parameter timeStepConstraint: (optional) set the method used to determine the integrator time step. "advective","oscillatory","min"
            arguments
                self WVModel {mustBeNonempty}
                options.deltaT (1,1) double {mustBePositive}
                options.cfl (1,1) double
                options.timeStepConstraint char {mustBeMember(options.timeStepConstraint,["advective","oscillatory","min"])} = "advective"
            end

            self.integratorOptions = options;

            self.createFixedTimeStepIntegrator();
            self.didSetupIntegrator = 1;
        end

        function resetFixedTimeStepIntegrator(self)
            self.integratorOptions = [];
            self.rk4Integrator = [];
            self.stepsTaken=0;    % number of RK4 steps/increments that have been made
            self.nSteps=inf;      % total number of expected RK4 increments to reach the model time requested by the user
            self.stepsPerOutput = [];      % number of RK4 steps between each output
            self.firstOutputStep = [];
        end

        function createFixedTimeStepIntegrator(self)
            if length(self.outputGroups) > 1
                error('The fixed time step integrator only supports a single output interval; you can only have one output group.');
            end
            self.outputInterval = self.outputGroupWithName(self.defaultOutputGroupName).outputInterval;
            if ~isempty(self.outputInterval)
                effectiveOutputInterval = self.outputInterval;
            else
                effectiveOutputInterval = self.finalIntegrationTime-self.t;
            end

            if isfield(self.integratorOptions,"deltaT")
                if isfield(self.integratorOptions,"cfl")
                    warning('deltaT was already set, ignoring cfl')
                end
                deltaT = effectiveOutputInterval/ceil(effectiveOutputInterval/self.integratorOptions.deltaT);
                if self.integratorOptions.deltaT ~= deltaT
                    warning('deltaT changed from %f to %f to match the outputInterval.',self.integratorOptions.deltaT,deltaT);
                end
            else
                if ~isfield(self.integratorOptions,"cfl")
                    self.integratorOptions.cfl = 0.25;
                end
                [deltaT,advectiveDT,oscillatoryDT] = self.timeStepForCFL(self.integratorOptions.cfl,effectiveOutputInterval);
                if strcmp(self.integratorOptions.timeStepConstraint,"advective")
                    deltaT = advectiveDT;
                    fprintf('Using the advective dt')        
                elseif strcmp(self.integratorOptions.timeStepConstraint,"oscillatory")
                    deltaT = oscillatoryDT;
                    fprintf('Using the oscillatory dt')
                elseif strcmp(self.integratorOptions.timeStepConstraint,"min")
                    deltaT = min(oscillatoryDT,advectiveDT);
                    fprintf('Using the min dt')
                end
                fprintf(': %.2f s (%d steps per output)\n',deltaT,round(effectiveOutputInterval/deltaT));
            end

            % Now set the initial conditions and point the integrator to
            % the correct flux function
            Y0 = self.initialConditionsCellArray();
            if isempty(Y0{1})
                error('Nothing to do! You must have set to linear dynamics, without floats, drifters or tracers.');
            end
            self.rk4Integrator = WVArrayIntegrator(@(t,y0) self.fluxAtTimeCellArray(t,y0),Y0,deltaT,currentTime=self.t);

            self.stepsPerOutput = round(effectiveOutputInterval/self.rk4Integrator.stepSize);
            self.firstOutputStep = round((self.initialOutputTime-self.t)/self.rk4Integrator.stepSize);
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
            U = self.wvt.uvMax;
            dx = (3/2)*(self.wvt.x(2)-self.wvt.x(1));

            advectiveDT = cfl*dx/U;
            oscillatoryDT = cfl*period;
            % A cfl of 1/12 for oscillatoryDT might be necessary for good numerical precision when advecting particles.

            if self.wvt.isBarotropic ~= 1
                W = self.wvt.w;
                W = W(:,:,2:end);
                ratio = abs((W./((3/2)*shiftdim(diff(self.wvt.z),-2))));
                dZW = 1/max(ratio(:));
                verticalAdvectiveDT = cfl*dZW;
                advectiveDT = min(verticalAdvectiveDT,advectiveDT);
                fprintf('dX/U = %.1f s (%.1f min). dZ/W = %.1f s (%.1f min). The highest frequency resolved IGW has period of %.1f s (%.1f min).\n', dx/U,dx/U/60, dZW,dZW/60,period,period/60);
            else
                fprintf('dX/U = %.1f s (%.1f min). The highest frequency resolved IGW has period of %.1f s (%.1f min).\n', dx/U,dx/U/60,period,period/60);
            end

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

        function integrateToTimeWithFixedTimeStep(self,finalTime)
            % Time step the model forward to the requested time.
            % - Topic: Integration
            arguments
                self WVModel {mustBeNonempty}
                finalTime (1,1) double
            end
            self.createFixedTimeStepIntegrator();

            self.finalIntegrationTime = finalTime;
            self.nSteps = ceil((self.finalIntegrationTime-self.t)/self.rk4Integrator.stepSize);               

            self.showIntegrationStartDiagnostics(self.finalIntegrationTime);
            while(self.t < self.finalIntegrationTime)
                self.integrateOneTimeStep();
                if self.didBlowUp == 1
                    return;
                end
                self.showIntegrationTimeDiagnostics(self.finalIntegrationTime);
            end
            self.showIntegrationFinishDiagnostics();
            self.finalIntegrationTime = [];
        end

        function Y0 = initialConditionsCellArray(self)
            Y0 = cell(1,1);
            n = 0;
            if self.linearDynamics == 0
                if self.wvt.hasWaveComponent == true
                    n=n+1;Y0{n} = self.wvt.Ap;
                    n=n+1;Y0{n} = self.wvt.Am;
                end
                if self.wvt.hasPVComponent == true
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


        function F = fluxAtTimeCellArray(self,t,y0)
            self.updateIntegratorValuesFromCellArray(t,y0);

            F = cell(1,1);
            n = 0;
            if self.linearDynamics == 0
                nlF = self.wvt.nFluxedComponents;
                [nlF{:}] = self.nonlinearFluxOperation.compute(self.wvt);
                if self.wvt.hasWaveComponent == true
                    n=n+1; F{n} = nlF{n};
                    n=n+1; F{n} = nlF{n};
                end
                if self.wvt.hasPVComponent == true
                    n=n+1; F{n} = nlF{n};
                end
            else

            end

            for iParticles=1:length(self.particle)
                p = self.particle{iParticles};
                if self.particle{iParticles}.fluxOp.isXYOnly
                    [F{n+1},F{n+2}] = self.particle{iParticles}.fluxOp.compute(self.wvt,p.x,p.y,p.z);
                    n=n+2;
                else
                    [F{n+1},F{n+2},F{n+3}] = self.particle{iParticles}.fluxOp.compute(self.wvt,p.x,p.y,p.z);
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

        function updateIntegratorValuesFromCellArray(self,t,y0)
            n=0;
            self.wvt.t = t;
            if self.linearDynamics == 0
                if self.wvt.hasWaveComponent == true
                    n=n+1; self.wvt.Ap = y0{n};
                    n=n+1; self.wvt.Am = y0{n};
                end
                if self.wvt.hasPVComponent == true
                    n=n+1; self.wvt.A0 = y0{n};
                end
            end

            for iParticles=1:length(self.particle)
                n=n+1; self.particle{iParticles}.x = y0{n};
                n=n+1; self.particle{iParticles}.y = y0{n};
                if ~self.particle{iParticles}.fluxOp.isXYOnly
                    n=n+1; self.particle{iParticles}.z = y0{n};
                end
            end

            for iTracer=1:length(self.tracerArray)
                n=n+1; self.tracerArray{iTracer} = self.y0{n};
            end
        end

        function modelTime = integrateOneTimeStep(self)
            % Ask the integrator to take one step forward, then record the
            % results.
            self.rk4Integrator.IncrementForward();

            % Rather than use the integrator time, which add floating point
            % numbers each time step, we multiple the steps taken by the
            % step size. This reduces rounding errors.
            self.stepsTaken = self.stepsTaken + 1;
            modelTime = self.initialTime + self.stepsTaken * self.rk4Integrator.stepSize;
            self.wvt.t = modelTime;

            self.updateIntegratorValuesFromCellArray(modelTime,self.rk4Integrator.currentY);

            if mod(self.stepsTaken - self.firstOutputStep,self.stepsPerOutput) == 0
                self.writeTimeStepToNetCDFFile();
            end

        end

        function modelTime = integrateToNextOutputTime(self)
            % Time step the model forward to the next output time
            % - Topic: Integration
            if isempty(self.outputInterval)
                fprintf('You did not set an output interval, so how could I integrateToNextOutputTime?\n');
                return;
            end

            self.openNetCDFFileForTimeStepping();
            if mod(self.stepsTaken - self.firstOutputStep,self.stepsPerOutput) == 0
                modelTime = self.integrateOneTimeStep;
                if self.didBlowUp == 1
                    self.recordNetCDFFileHistory(didBlowUp=1);
                    return;
                end
            end

            while( mod(self.stepsTaken - self.firstOutputStep,self.stepsPerOutput) ~= 0 )
                modelTime = self.integrateOneTimeStep;
                if self.didBlowUp == 1
                    self.recordNetCDFFileHistory(didBlowUp=1);
                    return;
                end
            end

            self.recordNetCDFFileHistory();
        end

    end
end