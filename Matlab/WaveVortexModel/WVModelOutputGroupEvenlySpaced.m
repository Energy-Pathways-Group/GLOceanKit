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

        initialTime
        finalTime
    end

    methods
        function self = WVModelOutputGroupEvenlySpaced(model,options)
            arguments
                model WVModel
                options.name {mustBeText}
                options.outputInterval (1,1) double {mustBePositive}
                options.initialTime (1,1) double = -Inf
                options.finalTime (1,1) double = Inf
            end
            self@WVModelOutputGroup(model,options.name);
            if ~isfield(options,"outputInterval")
                error("You must specify an output interval");
            end
            self.outputInterval = options.outputInterval;
            if options.initialTime == -Inf
                options.initialTime = model.wvt.t;
            end
            self.initialTime = options.initialTime;
            self.finalTime = options.finalTime;
        end

        function t = outputTimesForIntegrationPeriod(self,initialTime,finalTime)
            arguments (Input)
                self WVModelOutputGroup
                initialTime (1,1) double
                finalTime (1,1) double
            end
            arguments (Output)
                t (:,1) double
            end
            % Two possibilities here either:
            % 1) nothing is initialize, or 2) we already wrote to file
            if self.timeOfLastIncrementWrittenToFile == -Inf
                t = (self.initialTime:self.outputInterval:finalTime).';
            else
                t = ((self.timeOfLastIncrementWrittenToFile+self.outputInterval):self.outputInterval:finalTime).';
            end
            t(t<initialTime) = [];
            t(t<self.initialTime) = [];
            t(t>self.finalTime) = [];
        end

    end
end