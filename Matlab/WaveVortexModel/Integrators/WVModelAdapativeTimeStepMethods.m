classdef WVModelAdapativeTimeStepMethods < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    % properties (Abstract,GetAccess=public, SetAccess=public)
    %     nonlinearFluxOperation WVNonlinearFluxOperation
    % end
    properties (Abstract,GetAccess=public, SetAccess=protected)
        wvt
    end
    properties (Abstract) %(Access = protected)
        particle
        tracerArray
        linearDynamics
        didSetupIntegrator
        finalIntegrationTime
    end
    methods (Abstract)
        flag = didBlowUp(self)
        showIntegrationStartDiagnostics(self,finalTime)
        showIntegrationTimeDiagnostics(self,finalTime)
        showIntegrationFinishDiagnostics(self)
    end
    properties %(Access = protected)
        arrayLength
        arrayStartIndex
        arrayEndIndex
        odeOptions
        odeIntegrator

        % list of remaining times that need to be output to file
        outputTimes = [];
    end

    methods

        function setupAdaptiveTimeStepIntegrator(self,options)
            arguments
                self WVModel {mustBeNonempty}
                options.shouldShowIntegrationStats double {mustBeMember(options.shouldShowIntegrationStats,[0 1])} = 0
                options.integrator = @ode78
                options.absTolerance = 1e-6
                options.relTolerance = 1e-3;
                options.shouldUseScaledTolerance = 1;
                options.absToleranceA0 = 1e-10
                options.absToleranceApm = 1e-6
                options.absToleranceXY = 1e-1; % 100 km * 10^{-6}
                options.absToleranceZ = 1e-2;  
            end

            nArray = self.lengthOfFluxComponents;
            startIndex = 1;
            for i=1:length(nArray)
                self.arrayStartIndex(i) = startIndex;
                self.arrayEndIndex(i) = startIndex+nArray(i)-1;
                startIndex = self.arrayEndIndex(i)+1;
            end
            self.arrayLength = sum(nArray);

            function absTol = absoluteErrorTolerance()
                absTol = zeros(self.arrayLength,1);

                n = 0;
                if self.linearDynamics == 0
                    if self.wvt.hasWaveComponent == true
                        n=n+1;absTol(self.arrayStartIndex(n):self.arrayEndIndex(n)) = options.absToleranceApm;
                        n=n+1;absTol(self.arrayStartIndex(n):self.arrayEndIndex(n)) = options.absToleranceApm;
                    end
                    if self.wvt.hasPVComponent == true
                        n=n+1;absTol(self.arrayStartIndex(n):self.arrayEndIndex(n)) = options.absToleranceA0;
                    end
                end

                for iParticles=1:length(self.particle)
                    n=n+1;absTol(self.arrayStartIndex(n):self.arrayEndIndex(n)) = options.absToleranceXY;
                    n=n+1;absTol(self.arrayStartIndex(n):self.arrayEndIndex(n)) = options.absToleranceXY;
                    if ~self.particle{iParticles}.fluxOp.isXYOnly
                        n=n+1;absTol(self.arrayStartIndex(n):self.arrayEndIndex(n)) = options.absToleranceZ;
                    end
                end

                for iTracer=1:length(self.tracerArray)
                    n=n+1;absTol(self.arrayStartIndex(n):self.arrayEndIndex(n)) = 1e-5;
                end
            end

            function absTol = absoluteErrorToleranceAdaptive()
                absTol = zeros(self.arrayLength,1);

                alpha0 = options.absTolerance*sqrt(1./self.wvt.A0_TE_factor);
                alphapm = options.absTolerance*sqrt(1./self.wvt.Apm_TE_factor);
                alpha0(isinf(alpha0)) = 1;
                alphapm(isinf(alphapm)) = 1;

                n = 0;
                if self.linearDynamics == 0
                    if self.wvt.hasWaveComponent == true
                        n=n+1;absTol(self.arrayStartIndex(n):self.arrayEndIndex(n)) = alphapm;
                        n=n+1;absTol(self.arrayStartIndex(n):self.arrayEndIndex(n)) = alphapm;
                    end
                    if self.wvt.hasPVComponent == true
                        n=n+1;absTol(self.arrayStartIndex(n):self.arrayEndIndex(n)) = alpha0;
                    end
                end

                for iParticles=1:length(self.particle)
                    n=n+1;absTol(self.arrayStartIndex(n):self.arrayEndIndex(n)) = options.absToleranceXY;
                    n=n+1;absTol(self.arrayStartIndex(n):self.arrayEndIndex(n)) = options.absToleranceXY;
                    if ~self.particle{iParticles}.fluxOp.isXYOnly
                        n=n+1;absTol(self.arrayStartIndex(n):self.arrayEndIndex(n)) = options.absToleranceZ;
                    end
                end

                for iTracer=1:length(self.tracerArray)
                    n=n+1;absTol(self.arrayStartIndex(n):self.arrayEndIndex(n)) = 1e-5;
                end
            end

            self.odeOptions = odeset('OutputFcn',@self.timeStepIncrement);
            self.odeOptions = odeset(self.odeOptions,'RelTol',options.relTolerance);
            if options.shouldUseScaledTolerance == 1
                self.odeOptions = odeset(self.odeOptions,'AbsTol',absoluteErrorToleranceAdaptive);
            else
                self.odeOptions = odeset(self.odeOptions,'AbsTol',absoluteErrorTolerance);
            end
            self.odeOptions = odeset(self.odeOptions,'Refine',1); % must be set to 1
            if options.shouldShowIntegrationStats == 1
                self.odeOptions = odeset(self.odeOptions,'Stats','on');
            else
                self.odeOptions = odeset(self.odeOptions,'Stats','off');
            end

            self.odeIntegrator = options.integrator;

            self.didSetupIntegrator = 1;
        end

        function resetAdapativeTimeStepIntegrator(self)
            self.arrayLength = [];
            self.arrayStartIndex = [];
            self.arrayEndIndex = [];
            self.odeOptions = [];
            self.odeIntegrator = [];
        end

        function integrateToTimeWithAdaptiveTimeStep(self,finalTime)
            % Time step the model forward to the requested time.
            % - Topic: Integration
            arguments
                self WVModel {mustBeNonempty}
                finalTime (1,1) double
            end
            if ~isempty(self.ncfile)
                self.outputTimes = ((self.timeOfLastIncrementWrittenToFile+self.outputInterval):self.outputInterval:finalTime).';
                integratorTimes = self.outputTimes;
                if integratorTimes(1) ~= self.t
                    integratorTimes = cat(1,self.t,integratorTimes);
                end
                if integratorTimes(end) ~= finalTime
                    integratorTimes = cat(1,integratorTimes,finalTime);
                end
            else
                integratorTimes = [self.t finalTime];
            end

            self.finalIntegrationTime = finalTime;
            self.odeIntegrator(@(t,y) self.fluxArray(t,y),integratorTimes,self.initialConditionsArray,self.odeOptions);
            self.finalIntegrationTime = [];
        end

        function status = timeStepIncrement(self,t,y,flag)
            % Important notes:
            % because we set odeset(options,'Refine',1), t should only have
            % 1 value, other than for init and done. We are depending on
            % this behavior in the logic below.
            if strcmp(flag,'init')
                self.showIntegrationStartDiagnostics(self.finalIntegrationTime);
            elseif strcmp(flag,'done')
                self.showIntegrationFinishDiagnostics();
            else
                for iTime=1:length(t)
                        self.updateIntegratorValues(t(iTime),y(:,iTime))

                    if ~isempty(self.outputTimes) && abs(t(iTime) - self.outputTimes(1)) < eps
                        self.writeTimeStepToNetCDFFile();
                        self.outputTimes(1) = [];
                    end
                end
                self.showIntegrationTimeDiagnostics(self.finalIntegrationTime);
            end

            if self.didBlowUp == 1
                status = 1;
            else
                status = 0;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Integration: initial conditions, flux, and one time step
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function nArray = lengthOfFluxComponents(self)
            n = 0;
            if self.wvt.hasWaveComponent == true
                n=n+1; nArray(n) = numel(self.wvt.Ap);
                n=n+1; nArray(n) = numel(self.wvt.Am);
            end
            if self.wvt.hasPVComponent == true
                n=n+1; nArray(n) = numel(self.wvt.A0);
            end

            for iParticles=1:length(self.particle)
                p = self.particle{iParticles};
                n=n+1;nArray(n) = numel(p.x);
                n=n+1;nArray(n) = numel(p.y);
                if ~self.particle{iParticles}.fluxOp.isXYOnly
                    n=n+1;nArray(n) = numel(p.z);
                end
            end

            for i=1:length(self.tracerArray)
                n=n+1;nArray(n) = numel(self.tracerArray{i});
            end
        end

        function Y0 = initialConditionsArray(self)
            Y0 = zeros(self.arrayLength,1);

            n = 0;
            if self.linearDynamics == 0
                if self.wvt.hasWaveComponent == true
                    n=n+1;Y0(self.arrayStartIndex(n):self.arrayEndIndex(n)) = self.wvt.Ap(:);
                    n=n+1;Y0(self.arrayStartIndex(n):self.arrayEndIndex(n)) = self.wvt.Am(:);
                end
                if self.wvt.hasPVComponent == true
                    n=n+1;Y0(self.arrayStartIndex(n):self.arrayEndIndex(n)) = self.wvt.A0(:);
                end
            end

            for iParticles=1:length(self.particle)
                p = self.particle{iParticles};
                n=n+1;Y0(self.arrayStartIndex(n):self.arrayEndIndex(n)) = p.x(:);
                n=n+1;Y0(self.arrayStartIndex(n):self.arrayEndIndex(n)) = p.y(:);
                if ~self.particle{iParticles}.fluxOp.isXYOnly
                    n=n+1;Y0(self.arrayStartIndex(n):self.arrayEndIndex(n)) = p.z(:);
                end
            end

            for i=1:length(self.tracerArray)
                n=n+1;Y0(self.arrayStartIndex(n):self.arrayEndIndex(n)) = self.tracerArray{i}(:);
            end
        end

        function F = fluxArray(self,t,y0)
            self.updateIntegratorValues(t,y0)

            F = zeros(self.arrayLength,1);
            n = 0;
            if self.linearDynamics == 0
                nlF = cell(1,self.wvt.nFluxedComponents);
                [nlF{:}] = self.wvt.nonlinearFlux();
                if self.wvt.hasWaveComponent == true
                    n=n+1; F(self.arrayStartIndex(n):self.arrayEndIndex(n)) = nlF{n};
                    n=n+1; F(self.arrayStartIndex(n):self.arrayEndIndex(n)) = nlF{n};
                end
                if self.wvt.hasPVComponent == true
                    n=n+1; F(self.arrayStartIndex(n):self.arrayEndIndex(n)) = nlF{n};
                end
            else

            end

            for iParticles=1:length(self.particle)
                p = self.particle{iParticles};
                if self.particle{iParticles}.fluxOp.isXYOnly
                    [F(self.arrayStartIndex(n+1):self.arrayEndIndex(n+1)),F(self.arrayStartIndex(n+2):self.arrayEndIndex(n+2))] = self.particle{iParticles}.fluxOp.compute(self.wvt,p.x,p.y,p.z);
                    n=n+2;
                else
                    [F(self.arrayStartIndex(n+1):self.arrayEndIndex(n+1)),F(self.arrayStartIndex(n+2):self.arrayEndIndex(n+2)),F(self.arrayStartIndex(n+3):self.arrayEndIndex(n+3))] = self.particle{iParticles}.fluxOp.compute(self.wvt,p.x,p.y,p.z);
                    n=n+3;
                end
            end

            if ~isempty(self.tracerArray)
                for i=1:length(self.tracerArray)
                    phibar = self.wvt.transformFromSpatialDomainWithF(y0{n+1});
                    [~,Phi_x,Phi_y,Phi_z] = self.wvt.transformToSpatialDomainWithFAllDerivatives(phibar);
                    n=n+1;F(self.arrayStartIndex(n):self.arrayEndIndex(n)) = -self.wvt.u .* Phi_x - self.wvt.v.*Phi_y - self.wvt.w.*Phi_z;
                end
            end
        end


        function updateIntegratorValues(self,t,y0)
            n=0;
            self.wvt.t = t;
            if self.linearDynamics == 0
                if self.wvt.hasWaveComponent == true
                    n=n+1; self.wvt.Ap(:) = y0(self.arrayStartIndex(n):self.arrayEndIndex(n));
                    n=n+1; self.wvt.Am(:) = y0(self.arrayStartIndex(n):self.arrayEndIndex(n));
                end
                if self.wvt.hasPVComponent == true
                    n=n+1; self.wvt.A0(:) = y0(self.arrayStartIndex(n):self.arrayEndIndex(n));
                end
            end

            for iParticles=1:length(self.particle)
                n=n+1; self.particle{iParticles}.x(:) = y0(self.arrayStartIndex(n):self.arrayEndIndex(n));
                n=n+1; self.particle{iParticles}.y(:) = y0(self.arrayStartIndex(n):self.arrayEndIndex(n));
                if ~self.particle{iParticles}.fluxOp.isXYOnly
                    n=n+1; self.particle{iParticles}.z(:) = y0(self.arrayStartIndex(n):self.arrayEndIndex(n));
                end
            end

            for iTracer=1:length(self.tracerArray)
                n=n+1; self.tracerArray{iTracer}(:) = y0(self.arrayStartIndex(n):self.arrayEndIndex(n));
            end
        end

    end
end