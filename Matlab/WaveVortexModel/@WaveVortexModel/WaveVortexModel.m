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
    % - OpenNetCDFFileForTimeStepping should report expected file size

    % AAGH. Getting myself twisted in knots over the right API.
    % while ( tool.integrateToTime(finalTime) )
    %   tool.incrementForward()
    %   tool.WriteToFile()
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
    % Option 2: SetupIntegrator(deltaT)
    % Option 3: SetupIntegrator(deltaT,outputInterval)
    % Option 4: SetupIntegrator(deltaT,outputInterval,finalTime) -- alt: SetupIntegratorForFixedIntegrationTime
    % return estimated time steps?
    %
    % Integrate (called repeatedly):
    % Option 1: modelTime = integrateOneTimeStep()
    % Option 2: modelTime = integrateToNextOutputTime()
    % Option 3: modelTime = integrateToTime(futureTime)
    %
    % ShowIntegrationTimeDiagnostics( ???? )
    %
    %
    % Setup NetCDF
    % CreateNetCDFFileForModelOutput
    % AppendToExisting

    properties
        wvt             % WaveVortexTransform
        t=0             % current model time (in seconds)
        initialTime=0
    
        linearDynamics = 0

        % Variables integrated by this integration tool
        xFloat, yFloat, zFloat
        xDrifter, yDrifter, zDrifter
        tracers, tracerNames



        IMA0, IMAp, IMAm    % InteractionMasks
        EMA0, EMAp, EMAm    % EnergyMasks
        shouldAntiAlias = 1;

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
        variablesToWriteToFile = {}

        incrementsWrittenToFile
    end

    properties (SetAccess = private)
        particleIndexWithName % can't use a Map() because order matters
        particleFlux = {}
        particlePosition = {}

        tracerIndexWithName
        tracer = {}

        didSetupIntegrator=0
        initialConditionOnlyVariables
        timeSeriesVariables % internal use only
    end

    methods

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = WaveVortexModel(varargin)
            if isa(varargin{1},'WaveVortexTransform')
                WaveVortexTransform = varargin{1};
                if nargin > 1 && isa(varargin{2},"double")
                    self.t = varargin{2};
                else
                    self.t = 0;
                end
            elseif isa(varargin{1},'char' )
                existingModelOutput = varargin{1};

                restartIndex = Inf;
                shouldDoubleResolution = 0;
                if nargin > 1
                    restartIndex = varargin{2};
                end
                if nargin > 2
                    restartIndex = varargin{3};
                end

                nctool = WaveVortexModelNetCDFFile(existingModelOutput,'timeIndex',restartIndex);

                self.t = nctool.t;
                
                if (shouldDoubleResolution == 0)
                    WaveVortexTransform = nctool.wvm;
                else
                    WaveVortexTransform = nctool.wvm.WaveVortexTransformWithResolution(2*[nctool.wvm.Nx,nctool.wvm.Ny,nctool.wvm.nModes]);
                end

                % if there's existing model output, use that output interval
                time = ncread(existingModelOutput,'t');
                if length(time)>1
                    self.outputInterval = time(2)-time(1);
                end
            end

            self.initialOutputTime = self.t;
            self.initialTime = self.t;
            self.wvt = WaveVortexTransform;  

            % Allow all nonlinear interactions
            self.IMA0 = ones(self.wvt.Nk,self.wvt.Nl,self.wvt.nModes);
            self.IMAp = ones(self.wvt.Nk,self.wvt.Nl,self.wvt.nModes);
            self.IMAm = ones(self.wvt.Nk,self.wvt.Nl,self.wvt.nModes);

            % Allow energy fluxes at all modes
            self.EMA0 = ones(self.wvt.Nk,self.wvt.Nl,self.wvt.nModes);
            self.EMAp = ones(self.wvt.Nk,self.wvt.Nl,self.wvt.nModes);
            self.EMAm = ones(self.wvt.Nk,self.wvt.Nl,self.wvt.nModes);

            if self.shouldAntiAlias == 1
                self.disallowNonlinearInteractionsWithAliasedModes();
                self.freezeEnergyOfAliasedModes();
            end

            fluxVar(1) = ModelVariable('Fp',{'k','l','j'},'m/s2', 'non-linear flux into Ap with interaction and energy flux masks applied');
            fluxVar(2) = ModelVariable('Fm',{'k','l','j'},'m/s2', 'non-linear flux into Am with interaction and energy flux masks applied');
            fluxVar(3) = ModelVariable('F0',{'k','l','j'},'m/s', 'non-linear flux into A0 with interaction and energy flux masks applied');
            self.wvt.addModelOperation(ModelOperation('NonlinearFlux',fluxVar,@(wvt) self.NonlinearFluxWithMasks(wvt)));

            self.particleIndexWithName = containers.Map();
            self.tracerIndexWithName = containers.Map();
        end
        
        function SetNonlinearFlux(self,modelOp)
            arguments
                self WaveVortexModel {mustBeNonempty}
                modelOp ModelOperation {mustBeNonempty}
            end
            self.linearDynamics = 0;
        end

        function AddParticles(self,name,fluxOp,x,y,z)
            arguments
                self WaveVortexModel {mustBeNonempty}
                name char {mustBeNonempty}
                fluxOp ParticleFluxOperation {mustBeNonempty}
                x (:,1) double
                y (:,1) double
                z (:,1) double
            end
            
            n = length(self.particleIndexWithName) + 1;
            self.particleIndexWithName(name) = n;
            self.particleFlux{n} = fluxOp;
            self.particlePosition{n} = cat(2,x,y,z);
        end

        function [x,y,z] = ParticlePositions(self,name)
            p = self.particlePosition(self.particleIndexWithName(name));
            x = p(:,1);
            y = p(:,2);
            z = p(:,3);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Floats and drifters and tracer!
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function SetFloatPositions(self,x,y,z)
            floatFlux = ParticleFluxOperation('floatFlux',@(wvt,x,y,z) wvt.InterpolatedFieldAtPosition(x,y,z,'spline',wvt.u,wvt.v,wvt.w));
            self.AddParticles('float',floatFlux,x,y,z);
        end

        function [x,y,z] = FloatPositions(self)
            [x,y,z] = self.ParticlePosition('float');
        end

        function SetDrifterPositions(self,x,y,z)
            drifterFlux = ParticleFluxOperation('floatFlux',@(wvt,x,y,z) wvt.InterpolatedFieldAtPosition(x,y,z,'spline',wvt.u,wvt.v),0);
            self.AddParticles('drifter',drifterFlux,x,y,z);
        end

        function [x,y,z] = DrifterPositions(self)
            [x,y,z] = self.ParticlePosition('drifter');
        end

        function AddTracer(self,phi,name)
            n = length(self.tracerIndexWithName) + 1;
            self.tracerIndexWithName(name) = n;
            self.tracer{n} = phi;
        end

        function phi = Tracer(self,name)
            phi = self.tracer(self.tracerIndexWithName(name));
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Reduced interaction models
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function allowNonlinearInteractionsWithModes(self,Ap,Am,A0)
            self.IMA0 = or(self.IMA0,A0);
            self.IMAm = or(self.IMAm,Am);
            self.IMAp = or(self.IMAp,Ap);
        end

        function allowNonlinearInteractionsWithConstituents(self,constituents)
            [ApmMask,A0Mask] = self.wvt.MasksForFlowContinuents(constituents);
            self.IMA0 = self.IMA0 | A0Mask;
            self.IMAm = self.IMAm | ApmMask;
            self.IMAp = self.IMAp | ApmMask;
        end

        function disallowNonlinearInteractionsWithConstituents(self,constituents)
            [ApmMask,A0Mask] = self.wvt.MasksForFlowContinuents(constituents);
            self.IMA0 = self.IMA0 & ~A0Mask;
            self.IMAm = self.IMAm & ~ApmMask;
            self.IMAp = self.IMAp & ~ApmMask;
        end

        function disallowNonlinearInteractionsWithAliasedModes(self)
            % Uses the 2/3 rule to prevent aliasing of Fourier modes.
            % The reality is that the vertical modes will still alias.
            % http://helper.ipam.ucla.edu/publications/mtws1/mtws1_12187.pdf
            self.shouldAntiAlias = 1;
            
            AntiAliasMask = self.wvt.MaskForAliasedModes();
            self.IMA0 = self.IMA0 & ~AntiAliasMask;
            self.IMAm = self.IMAm & ~AntiAliasMask;
            self.IMAp = self.IMAp & ~AntiAliasMask;
        end

        function unfreezeEnergyOfConstituents(self,constituents)
            [ApmMask,A0Mask] = self.wvt.MasksForFlowContinuents(constituents);
            self.EMA0 = self.EMA0 | A0Mask;
            self.EMAm = self.EMAm | ApmMask;
            self.EMAp = self.EMAp | ApmMask;
        end

        function freezeEnergyOfConstituents(self,constituents)
            [ApmMask,A0Mask] = self.wvt.MasksForFlowContinuents(constituents);
            self.EMA0 = self.EMA0 & ~A0Mask;
            self.EMAm = self.EMAm & ~ApmMask;
            self.EMAp = self.EMAp & ~ApmMask;
        end
        
        function freezeEnergyOfAliasedModes(self)
            % In addition to disallowing interaction to occur between modes
            % that are aliased, you may actually want to disallow energy to
            % even enter the aliased modes.
            AntiAliasMask = self.wvt.MaskForAliasedModes();
            self.EMA0 = self.EMA0 & ~AntiAliasMask;
            self.EMAm = self.EMAm & ~AntiAliasMask;
            self.EMAp = self.EMAp & ~AntiAliasMask;
        end

        function clearEnergyFromAliasedModes(self)
            % In addition to disallowing interaction to occur between modes
            % that are aliased, you may actually want to disallow energy to
            % even enter the aliased modes.
            AntiAliasMask = self.wvt.MaskForAliasedModes();
            self.wvt.A0 = self.wvt.A0 .* ~AntiAliasMask;
            self.wvt.Am = self.wvt.Am .* ~AntiAliasMask;
            self.wvt.Ap = self.wvt.Ap .* ~AntiAliasMask;
        end

        function flag = IsAntiAliased(self)
            AntiAliasMask = self.MaskForAliasedModes();

            % check if there are zeros at all the anti-alias indices
            flag = all((~self.IMA0 & AntiAliasMask) == AntiAliasMask,'all');
            flag = flag & all((~self.IMAm & AntiAliasMask) == AntiAliasMask,'all');
            flag = flag & all((~self.IMAp & AntiAliasMask) == AntiAliasMask,'all');
            flag = flag & all((~self.EMA0 & AntiAliasMask) == AntiAliasMask,'all');
            flag = flag & all((~self.EMAp & AntiAliasMask) == AntiAliasMask,'all');
            flag = flag & all((~self.EMAm & AntiAliasMask) == AntiAliasMask,'all');
        end

        
        function [omega,k,l] = addForcingWaveModes(self,kModes,lModes,jModes,phi,U,signs)
            [omega,k,l,kIndex,lIndex,jIndex,signIndex] = self.wvt.AddGriddedWavesWithWavemodes(kModes,lModes,jModes,phi,U,signs);
            for iMode=1:length(kModes)
                if (signIndex(iMode) == 1 || (kIndex(iMode) == 1 && lIndex(iMode) == 1) )
                    self.EMAp( kIndex(iMode),lIndex(iMode),jIndex(iMode)) = 0;
                    self.EMAp = WaveVortexTransform.MakeHermitian(self.EMAp);
                end

                if (signIndex(iMode) == -1 || (kIndex(iMode) == 1 && lIndex(iMode) == 1) )
                    self.EMAm( kIndex(iMode),lIndex(iMode),jIndex(iMode)) = 0;
                    self.EMAm = WaveVortexTransform.MakeHermitian(self.EMAm);
                end
            end
        end

        function stirWithConstituents(self,constituents)

        end

        function [Fp,Fm,F0] = NonlinearFluxWithMasks(self,wvt)
            fprintf('Custom');
            % Apply operator T_\omega---defined in (C2) in the manuscript

            % We also apply the interaction masks (IMA*) and energy masks
            % (EMA*). These *could* be precomputed multiplied directly into
            % the coefficients. A quick speed test shows this about a 5%
            % performance hit by not doing this.
            phase = exp(wvt.iOmega*(wvt.t-wvt.t0));
            Ap = self.IMAp .* wvt.Ap .* phase;
            Am = self.IMAm .* wvt.Am .* conj(phase);
            A0 = self.IMA0 .* wvt.A0;

            % Apply operator S---defined in (C4) in the manuscript
            Ubar = wvt.UAp.*Ap + wvt.UAm.*Am + wvt.UA0.*A0;
            Vbar = wvt.VAp.*Ap + wvt.VAm.*Am + wvt.VA0.*A0;
            Wbar = wvt.WAp.*Ap + wvt.WAm.*Am;
            Nbar = wvt.NAp.*Ap + wvt.NAm.*Am + wvt.NA0.*A0;

            % Finishing applying S, but also compute derivatives at the
            % same time
            [U,Ux,Uy,Uz] = wvt.TransformToSpatialDomainWithFAllDerivatives(Ubar);
            [V,Vx,Vy,Vz] = wvt.TransformToSpatialDomainWithFAllDerivatives(Vbar);
            W = wvt.TransformToSpatialDomainWithG(Wbar);
            [ETA,ETAx,ETAy,ETAz] = wvt.TransformToSpatialDomainWithGAllDerivatives(Nbar);

            % Compute the nonlinear terms in the spatial domain
            % (pseudospectral!)
            uNL = -U.*Ux - V.*Uy - W.*Uz;
            vNL = -U.*Vx - V.*Vy - W.*Vz;
            nNL = -U.*ETAx - V.*ETAy - W.*(ETAz + ETA.*shiftdim(wvt.dLnN2,-2));

            % Now apply the operator S^{-1} and then T_\omega^{-1}
            uNLbar = wvt.TransformFromSpatialDomainWithF(uNL);
            vNLbar = wvt.TransformFromSpatialDomainWithF(vNL);
            nNLbar = wvt.TransformFromSpatialDomainWithG(nNL);

            Fp = self.EMAp .* (wvt.ApU.*uNLbar + wvt.ApV.*vNLbar + wvt.ApN.*nNLbar) .* conj(phase);
            Fm = self.EMAm .* (wvt.AmU.*uNLbar + wvt.AmV.*vNLbar + wvt.AmN.*nNLbar) .* phase;
            F0 = self.EMA0 .* (wvt.A0U.*uNLbar + wvt.A0V.*vNLbar + wvt.A0N.*nNLbar);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Integration loop
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = IntegrateToTime(self,finalTime,cfl)
            if self.didSetupIntegrator ~= 1
                if nargin < 3 || isempty(cfl)
                    cfl=0.5;
                end
                deltaT = self.wvt.TimeStepForCFL(cfl,self.outputInterval);
                self.SetupIntegrator(deltaT,self.outputInterval,finalTime);
            end
   
            self.OpenNetCDFFileForTimeStepping();
            
            while(self.t < finalTime)
                
                self.integrateToNextOutputTime();

                self.WriteTimeStepToNetCDFFile();
            end

%             self.CloseNetCDFFile();
        end

        function varargout = SetupIntegrator(self,deltaT,outputInterval,finalTime)
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
                deltaT = self.wvt.TimeStepForCFL(0.5);
            elseif didSetDeltaT == 0 && didSetOutputInterval == 1
                deltaT = self.wvt.TimeStepForCFL(0.5,outputInterval);
            end
            
            % Now set the initial conditions and point the integrator to
            % the correct flux function
            Y0 = self.InitialConditionsArray();
            if isempty(Y0{1})
                error('Nothing to do! You must have set to linear dynamics, without floats, drifters or tracers.');
            end
            self.integrator = ArrayIntegrator(@(t,y0) self.FluxAtTime(t,y0),Y0,deltaT);
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
        
        function Y0 = InitialConditionsArray(self)
            Y0 = cell(1,1);
            n = 0;
            if self.linearDynamics == 0
                n=n+1;Y0{n} = self.wvt.Ap;
                n=n+1;Y0{n} = self.wvt.Am;
                n=n+1;Y0{n} = self.wvt.A0;
            end

            if ~isempty(self.xFloat)
                n=n+1;Y0{n} = self.xFloat;
                n=n+1;Y0{n} = self.yFloat;
                n=n+1;Y0{n} = self.zFloat;
            end

            if ~isempty(self.xDrifter)
                n=n+1;Y0{n} = self.xDrifter;
                n=n+1;Y0{n} = self.yDrifter;
            end

            if ~isempty(self.tracers)
                for i=1:length(self.tracers)
                    n=n+1;Y0{n} = self.tracers{i};
                end
            end
        end

        function modelTime = integrateOneTimeStep(self)
            %self.ShowIntegrationTimeDiagnostics(self.stepsTaken);

            self.integrator.IncrementForward();
            n=0;
            if self.linearDynamics == 0
                n=n+1; self.wvt.Ap = self.integrator.currentY{n};
                n=n+1; self.wvt.Am = self.integrator.currentY{n};
                n=n+1; self.wvt.A0 = self.integrator.currentY{n};
            end

            for iParticles=1:length(self.particleFlux)
                p = self.particlePosition(iParticles);
                if self.particleFlux(iParticles).xyOnly
                    
                else

                end
            end

            if ~isempty(self.xFloat)
                n=n+1; self.xFloat = self.integrator.currentY{n};
                n=n+1; self.yFloat = self.integrator.currentY{n};
                n=n+1; self.zFloat = self.integrator.currentY{n};
            end

            if ~isempty(self.xDrifter)
                n=n+1; self.xDrifter = self.integrator.currentY{n};
                n=n+1; self.yDrifter = self.integrator.currentY{n};
            end

            if ~isempty(self.tracers)
                for iTracer=1:length(self.tracers)
                    n=n+1; self.tracers{iTracer} = self.integrator.currentY{n};
                end
            end

            % Rather than use the integrator time, which add floating point
            % numbers each time step, we multiple the steps taken by the
            % step size. This reduces rounding errors.
            self.stepsTaken = self.stepsTaken + 1;
            modelTime = self.initialTime + self.stepsTaken * self.integrator.stepSize;
            self.t = modelTime;
            if mod(self.stepsTaken - self.firstOutputStep,self.stepsPerOutput) == 0
                self.outputIndex = self.outputIndex + 1;
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

        function F = FluxAtTime(self,t,y0)
            F = cell(1,1);
            n = 0;
            self.wvt.t = t;
            if self.linearDynamics == 0
                [Fp,Fm,F0] = self.wvt.NonlinearFlux;
                n=n+1;F{n} = Fp;
                n=n+1;F{n} = Fm;
                n=n+1;F{n} = F0;
            end

            if ~isempty(self.xFloat)
                [Fx,Fy,Fz] = self.wvt.InterpolatedFieldAtPosition(y0{n+1},y0{n+2},y0{n+3},'spline',U,V,W);
                n=n+1;F{n} = Fx;
                n=n+1;F{n} = Fy;
                n=n+1;F{n} = Fz;
            end

            if ~isempty(self.xDrifter)
                [Fx,Fy] = self.wvt.InterpolatedFieldAtPosition(y0{n+1},y0{n+2},self.zDrifter,'spline',U,V);
                n=n+1;F{n} = Fx;
                n=n+1;F{n} = Fy;
            end

            if ~isempty(self.tracers)
                for i=1:length(self.tracers)
                    phibar = self.wvt.TransformFromSpatialDomainWithF(y0{n+1});
                    [~,Phi_x,Phi_y,Phi_z] = self.wvt.TransformToSpatialDomainWithFAllDerivatives(phibar);
                    n=n+1;F{n} = -U.*Phi_x - V.*Phi_y - W.*Phi_z;
                end
            end
        end

        function ShowIntegrationTimeDiagnostics(self,integratorIncrement)
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

        function [deltaT,advectiveDT,oscillatoryDT] = TimeStepForCFL(self, cfl, outputInterval)
            % Return the time step (in seconds) to maintain the given cfl condition.
            % If the cfl condition is not given, 0.25 will be assumed.
            % If outputInterval is given, the time step will be rounded to evenly
            % divide the outputInterval.
            if nargin == 1
                cfl = 0.25;
            end

            omega = self.wvt.Omega;
            period = 2*pi/max(abs(omega(:)));
            [u,v] = self.wvt.VelocityFieldAtTime(0.0);
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

        function ncfile = CreateNetCDFFileForModelOutput(self,modelOutputFile,overwriteExisting)
            shouldOverwriteExisting = 0;
            if strcmp(overwriteExisting,'OVERWRITE_EXISTING')
                shouldOverwriteExisting = 1;
            end
            ncfile = WaveVortexModelNetCDFFile(self.wvt,modelOutputFile,'shouldOverwriteExisting',shouldOverwriteExisting);            
            self.ncfile = ncfile;
        end

        function OpenNetCDFFileForTimeStepping(self)
            if ~isempty(self.ncfile)
    
                % Gather up all the field variables and sort them into bins
                % on whether or not the user wants the entire time series
                % written, or just the initial conditions.
                self.timeSeriesVariables = {};
                self.initialConditionOnlyVariables = {};
                if self.shouldWriteA0 == ShouldWrite.timeSeries
                    self.timeSeriesVariables{end+1} = 'A0';
                elseif self.shouldWriteA0 == ShouldWrite.initialConditions
                    self.initialConditionOnlyVariables{end+1} = 'A0';
                end
                if self.shouldWriteAp == ShouldWrite.timeSeries
                    self.timeSeriesVariables{end+1} = 'Ap';
                elseif self.shouldWriteAp == ShouldWrite.initialConditions
                    self.initialConditionOnlyVariables{end+1} = 'Ap';
                end
                if self.shouldWriteAm == ShouldWrite.timeSeries
                    self.timeSeriesVariables{end+1} = 'Am';
                elseif self.shouldWriteAm == ShouldWrite.initialConditions
                    self.initialConditionOnlyVariables{end+1} = 'Am';
                end

                if self.shouldWriteU == ShouldWrite.timeSeries
                    self.timeSeriesVariables{end+1} = 'u';
                elseif self.shouldWriteU == ShouldWrite.initialConditions
                    self.initialConditionOnlyVariables{end+1} = 'u';
                end
                if self.shouldWriteV == ShouldWrite.timeSeries
                    self.timeSeriesVariables{end+1} = 'v';
                elseif self.shouldWriteV == ShouldWrite.initialConditions
                    self.initialConditionOnlyVariables{end+1} = 'v';
                end
                if self.shouldWriteW == ShouldWrite.timeSeries
                    self.timeSeriesVariables{end+1} = 'w';
                elseif self.shouldWriteW == ShouldWrite.initialConditions
                    self.initialConditionOnlyVariables{end+1} = 'w';
                end

                if self.shouldWriteEta == ShouldWrite.timeSeries
                    self.timeSeriesVariables{end+1} = 'eta';
                elseif self.shouldWriteEta == ShouldWrite.initialConditions
                    self.initialConditionOnlyVariables{end+1} = 'eta';
                end
                if self.shouldWriteP == ShouldWrite.timeSeries
                    self.timeSeriesVariables{end+1} = 'p';
                elseif self.shouldWriteP == ShouldWrite.initialConditions
                    self.initialConditionOnlyVariables{end+1} = 'p';
                end
                if self.shouldWriteRho == ShouldWrite.timeSeries
                    self.timeSeriesVariables{end+1} = 'rho_prime';
                elseif self.shouldWriteRho == ShouldWrite.initialConditions
                    self.initialConditionOnlyVariables{end+1} = 'rho_prime';
                end

                if ~isempty(self.timeSeriesVariables)
                    self.ncfile.InitializeStorageForVariableFieldsTimeSeries(self.timeSeriesVariables);
                end
                if ~isempty(self.initialConditionOnlyVariables)
                    self.ncfile.WriteInitialVariableFields(self.initialConditionOnlyVariables);
                end

                if ~isempty(self.xFloat)
                    if self.shouldWriteFloats == ShouldWrite.timeSeries
                        self.ncfile.InitializeFloatStorageForTimeSeries(length(self.xFloat));
                    elseif self.shouldWriteFloats == ShouldWrite.initialConditions
                        self.ncfile.WriteInitialFloatPositions(self.xFloat,self.yFloat,self.zFloat);
                    end
                end

                if ~isempty(self.xDrifter)
                    if self.shouldWriteDrifters == ShouldWrite.timeSeries
                        self.ncfile.InitializeDrifterStorageForTimeSeries(length(self.xDrifter));
                    elseif self.shouldWriteDrifters == ShouldWrite.initialConditions
                        self.ncfile.WriteInitialDrifterPositions(self.xDrifter,self.yDrifter,self.zDrifter);
                    end
                end

                if self.shouldWriteTracers == ShouldWrite.timeSeries
                    for iTracer = 1:length(self.tracers)
                        self.ncfile.InitializeTracerStorageWithName(self.tracerNames{iTracer});
                    end
                elseif self.shouldWriteTracers == ShouldWrite.initialConditions
                    for iTracer = 1:length(self.tracers)
                        self.ncfile.WriteInitialTracerWithName(self.tracerNames{iTracer},self.tracers{iTracer});
                    end
                end

                if self.shouldWriteFlowConstituentEnergetics == ShouldWrite.timeSeries
                    self.ncfile.InitializeEnergeticsStorageForTimeSeries();
                elseif self.shouldWriteFlowConstituentEnergetics == ShouldWrite.initialConditions
                    self.ncfile.WriteInitialEnergetics();
                end

                if self.shouldWriteFlowConstituentEnergetics2D == ShouldWrite.timeSeries
                    self.ncfile.InitializeEnergeticsKJStorageForTimeSeries();
                elseif self.shouldWriteFlowConstituentEnergetics2D == ShouldWrite.initialConditions
                    self.ncfile.WriteInitialEnergeticsKJ();
                end
            else
                if isempty(self.ncfile.ncid)
                    self.ncfile.open();
                end
            end

            self.incrementsWrittenToFile = 0;

%             outputTimes = self.t + self.outputIncrements*self.integrator.stepSize;

            % Save the initial conditions
            self.WriteTimeStepToNetCDFFile();         
        end


        function WriteTimeStepToNetCDFFile(self)
            if ( ~isempty(self.ncfile) && mod(self.stepsTaken - self.firstOutputStep,self.stepsPerOutput) == 0 )
                self.ncfile.WriteTimeAtIndex(self.outputIndex,self.t);

                self.ncfile.WriteVariableFieldsAtTimeIndex(self.outputIndex,self.timeSeriesVariables{:});

                if self.shouldWriteFloats == ShouldWrite.timeSeries
                    self.ncfile.WriteFloatPositionsAtTimeIndex(self.outputIndex,self.xFloat,self.yFloat,self.zFloat);
                end
                if self.shouldWriteDrifters == ShouldWrite.timeSeries
                    self.ncfile.WriteDrifterPositionsTimeAtIndex(self.outputIndex,self.xDrifter,self.yDrifter,self.zDrifter);
                end
                if self.shouldWriteTracers == ShouldWrite.timeSeries
                    for iTracer = 1:length(self.tracers)
                        self.ncfile.WriteTracerWithNameTimeAtIndex(self.outputIndex,self.tracerNames{iTracer},self.tracers{iTracer});
                    end
                end

                if self.shouldWriteFlowConstituentEnergetics == ShouldWrite.timeSeries
                    self.ncfile.WriteEnergeticsAtTimeIndex(self.outputIndex);
                end
                if self.shouldWriteFlowConstituentEnergetics2D == ShouldWrite.timeSeries
                    self.ncfile.WriteEnergeticsKJAtTimeIndex(self.outputIndex);
                end

                self.incrementsWrittenToFile = self.incrementsWrittenToFile + 1;
            end
        end

        function CloseNetCDFFile(self)
            if ~isempty(self.ncfile)
                fprintf('Ending simulation. Wrote %d time points to file\n',self.incrementsWrittenToFile);
                self.ncfile.close();
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Integration
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        

    end
end