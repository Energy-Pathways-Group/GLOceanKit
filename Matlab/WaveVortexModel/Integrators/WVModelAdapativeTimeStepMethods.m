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
    end

    methods

        function setupAdaptiveTimeStepIntegrator(self,options)
            arguments
                self WVModel {mustBeNonempty}
                options.shouldShowIntegrationStats logical = false
                options.integrator = @ode78
                options.relTolerance = 1e-3;
            end

            nArray = self.lengthOfFluxComponents;
            startIndex = 1;
            for i=1:length(nArray)
                self.arrayStartIndex(i) = startIndex;
                self.arrayEndIndex(i) = startIndex+nArray(i)-1;
                startIndex = self.arrayEndIndex(i)+1;
            end
            self.arrayLength = sum(nArray);

            self.odeOptions = odeset('OutputFcn',@self.timeStepIncrementArray);
            self.odeOptions = odeset(self.odeOptions,'InitialStep',self.timeStepForCFL(0.5));
            self.odeOptions = odeset(self.odeOptions,'RelTol',options.relTolerance);
            self.odeOptions = odeset(self.odeOptions,'AbsTol',self.absErrorToleranceArray);
            self.odeOptions = odeset(self.odeOptions,'Refine',1); % must be set to 1
            if options.shouldShowIntegrationStats
                self.odeOptions = odeset(self.odeOptions,'Stats','on');
            else
                self.odeOptions = odeset(self.odeOptions,'Stats','off');
            end

            self.odeIntegrator = options.integrator;
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
            
            % The function call here is stupid, because it is not obvious
            % that callign outputTimesForIntegrationPeriod actually has the
            % side-effect of setting up the run
            integratorTimes = self.outputTimesForIntegrationPeriod(self.t,finalTime);
            arrayfun( @(outputFile) outputFile.writeTimeStepToOutputFile(self.t), self.outputFiles);

            self.finalIntegrationTime = finalTime;
            self.odeIntegrator(@(t,y) self.fluxArray(t,y),integratorTimes,self.initialConditionsArray,self.odeOptions);
            self.finalIntegrationTime = [];
        end

        function status = timeStepIncrementArray(self,t,y,flag)
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
                for iTime=1:length(t)
                    self.updateIntegratorValuesFromArray(t(iTime),y(:,iTime))
                    self.writeTimeStepToNetCDFFile(t(iTime));
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
        % Conversion from cell array to linear array
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function nArray = lengthOfFluxComponents(self)
            nArray = [];
            for i = 1:length(self.fluxedObservingSystems)
                nArray = cat(1,nArray,self.fluxedObservingSystems(i).lengthOfFluxComponents);
            end
        end

        function Y0_array = absErrorToleranceArray(self)
            Y0_array = zeros(self.arrayLength,1);
            Y0_cell = self.absErrorToleranceCellArray;

            for n=1:self.nFluxComponents
                Y0_array(self.arrayStartIndex(n):self.arrayEndIndex(n)) = Y0_cell{n}(:);
            end
        end

        function Y0_array = initialConditionsArray(self)
            Y0_array = zeros(self.arrayLength,1);
            Y0_cell = self.initialConditionsCellArray;

            for n=1:self.nFluxComponents
                Y0_array(self.arrayStartIndex(n):self.arrayEndIndex(n)) = Y0_cell{n}(:);
            end
        end

        function F_array = fluxArray(self,t,Y0_array)
            F_array = zeros(self.arrayLength,1);
            Y0_cell = cell(self.nFluxComponents,1);
            for n=1:self.nFluxComponents
                Y0_cell{n} = Y0_array(self.arrayStartIndex(n):self.arrayEndIndex(n));
            end

            F_cell = self.fluxAtTimeCellArray(t,Y0_cell);

            for n=1:self.nFluxComponents
                F_array(self.arrayStartIndex(n):self.arrayEndIndex(n)) = F_cell{n}(:);
            end
        end

        function updateIntegratorValuesFromArray(self,t,Y0_array)
            Y0_cell = cell(self.nFluxComponents,1);
            for n=1:self.nFluxComponents
                Y0_cell{n} = Y0_array(self.arrayStartIndex(n):self.arrayEndIndex(n));
            end
            self.updateIntegratorValuesFromCellArray(t,Y0_cell);
        end
    end
end