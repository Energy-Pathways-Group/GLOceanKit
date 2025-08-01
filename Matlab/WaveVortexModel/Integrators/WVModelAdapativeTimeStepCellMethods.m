classdef WVModelAdapativeTimeStepCellMethods < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    % properties (Abstract,GetAccess=public, SetAccess=public)
    %     nonlinearFluxOperation WVNonlinearFluxOperation
    % end
    properties (Abstract,GetAccess=public, SetAccess=protected)
        wvt
    end
    properties (Abstract) %(Access = protected)
        finalIntegrationTime
        odeOptions
        odeIntegrator
    end
    methods (Abstract)
        flag = didBlowUp(self)
        showIntegrationStartDiagnostics(self,finalTime)
        showIntegrationTimeDiagnostics(self,finalTime)
        showIntegrationFinishDiagnostics(self)

    end

    methods

        function setupAdaptiveTimeStepCellIntegrator(self,options)
            arguments
                self WVModel {mustBeNonempty}
                options.shouldShowIntegrationStats logical = false
                options.integrator = @ode45_cell
                options.relTolerance = 1e-3;
            end

            self.odeOptions = odeset('OutputFcn',@self.timeStepIncrement);
            self.odeOptions = odeset(self.odeOptions,'InitialStep',self.timeStepForCFL(0.5));
            self.odeOptions = odeset(self.odeOptions,'RelTol',options.relTolerance);
            self.odeOptions = odeset(self.odeOptions,'AbsTol',self.absErrorToleranceCellArray);
            self.odeOptions = odeset(self.odeOptions,'Refine',1); % must be set to 1
            if options.shouldShowIntegrationStats
                self.odeOptions = odeset(self.odeOptions,'Stats','on');
            else
                self.odeOptions = odeset(self.odeOptions,'Stats','off');
            end

            self.odeIntegrator = options.integrator;
        end

        function integrateToTimeWithAdaptiveTimeStepCell(self,finalTime)
            % Time step the model forward to the requested time.
            % - Topic: Integration
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
            self.odeIntegrator(@(t,y) self.fluxAtTimeCellArray(t,y),integratorTimes,self.initialConditionsCellArray,self.odeOptions);
            self.finalIntegrationTime = [];
        end
    end
end