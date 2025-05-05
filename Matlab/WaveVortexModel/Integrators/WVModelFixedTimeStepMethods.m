classdef WVModelFixedTimeStepMethods < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    % properties (Abstract,GetAccess=public, SetAccess=public)
    %     nonlinearFluxOperation
    % end
    properties (Abstract,GetAccess=public, SetAccess=protected)
        wvt
    end

    properties
        integratorOptions
        rk4Integrator
    end
    properties (Abstract) %(Access = protected)
        finalIntegrationTime
    end
    methods (Abstract)
        flag = didBlowUp(self)
        showIntegrationStartDiagnostics(self,finalTime)
        showIntegrationTimeDiagnostics(self,finalTime)
        showIntegrationFinishDiagnostics(self)
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
                options.timeStepConstraint char {mustBeMember(options.timeStepConstraint,["advective","oscillatory","min"])} = "min"
            end

            self.integratorOptions = options;
        end

        function resetFixedTimeStepIntegrator(self)
            self.integratorOptions = [];
            self.rk4Integrator = [];
        end

        function [deltaT,advectiveDT,oscillatoryDT] = timeStepForCFL(self, cfl, outputInterval)
            % Return the time step (in seconds) to maintain the given cfl condition.
            % If the cfl condition is not given, 0.25 will be assumed.
            % If outputInterval is given, the time step will be rounded to evenly
            % divide the outputInterval.
            if nargin == 1
                cfl = 0.25;
            end

            U = self.wvt.uvMax;
            dx = self.wvt.effectiveHorizontalGridResolution;
            advectiveDT = cfl*dx/U;
            oscillatoryDT = Inf;
            period = Inf;
            if self.wvt.hasWaveComponent == true
                omega = self.wvt.Omega;
                period = 2*pi/max(abs(omega(:)));
                oscillatoryDT = cfl*period;
            end
            % A cfl of 1/12 for oscillatoryDT might be necessary for good numerical precision when advecting particles.

            if self.wvt.hasVariableWithName("w")
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
            if isfield(self.integratorOptions,"deltaT")
                deltaT = self.integratorOptions.deltaT;
            else
                if ~isfield(self.integratorOptions,"cfl")
                    self.integratorOptions.cfl = 0.25;
                end
                [deltaT,advectiveDT,oscillatoryDT] = self.timeStepForCFL(self.integratorOptions.cfl);
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
            end

            % The function call here is stupid, because it is not obvious
            % that callign outputTimesForIntegrationPeriod actually has the
            % side-effect of setting up the run
            integratorTimes = self.outputTimesForIntegrationPeriod(self.t,finalTime);
            arrayfun( @(outputFile) outputFile.writeTimeStepToOutputFile(self.t), self.outputFiles);

            self.finalIntegrationTime = finalTime;
            Y0 = self.initialConditionsCellArray();
            self.rk4Integrator = WVArrayIntegrator(@(t,y0) self.fluxAtTimeCellArray(t,y0),integratorTimes,Y0,deltaT,OutputFcn=@self.timeStepIncrement);
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
                % fprintf("nFlux computations: " + string(self.nFluxComputations) + "\n");
            else
                % these are the 'interpolation' times, between the actual
                % time step
                self.updateIntegratorValuesFromCellArray(t,y)
                self.writeTimeStepToNetCDFFile(t);
                self.showIntegrationTimeDiagnostics(self.finalIntegrationTime);
            end

            if self.didBlowUp == 1
                status = 1;
            else
                status = 0;
            end
        end

    end
    
end