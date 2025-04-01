classdef WVModelOutputGroupEvenlySpaced < WVModelOutputGroup
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        % Model output interval (seconds)
        % - Topic: Integration
        % This property is optionally set when calling setupIntegrator. If
        % set, it will allow you to call -integrateToNextOutputTime and, if
        % a NetCDF file is set for output, it will set the interval at
        % which time steps are written to file.
        outputInterval = [] % (1,1) double

        t0 = 0;
    end

    methods
        function self = WVModelOutputGroupEvenlySpaced(model,options)
            arguments
                model WVModel
                options.name {mustBeText}
                options.outputInterval (1,1) double {mustBePositive}
            end
            self@WVModelOutputGroup(model,options.name);
            self.outputInterval = options.outputInterval;
        end

        function t = outputTimesForIntegrationPeriod(self,initialTime,finalTime)
            if self.timeOfLastIncrementWrittenToFile == -Inf
                error('self.timeOfLastIncrementWrittenToFile not initialized');
            end
            t = ((self.timeOfLastIncrementWrittenToFile+self.outputInterval):self.outputInterval:finalTime).';
            t(t<initialTime) = [];
        end

    end
end