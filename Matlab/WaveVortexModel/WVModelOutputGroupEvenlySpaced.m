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
            self@WVModelOutputGroup(model,name=options.name);
            if ~isfield(options,"name")
                error("You must specify an output group name");
            end
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

        function aString = description(self)
            if self.finalTime == Inf
                aString = "evenly spaced output with an interval of " + string(self.outputInterval) + " starting at time t=" + string(self.initialTime) + " and continuing indefinitely";
            else
                aString = "evenly spaced output with an interval of " + string(self.outputInterval) + " start at time t=" + string(self.initialTime) + "and ending at time t=" + string(self.finalTime);
            end
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
            if self.timeOfLastIncrementWrittenToGroup == -Inf
                t = (self.initialTime:self.outputInterval:finalTime).';
            else
                t = ((self.timeOfLastIncrementWrittenToGroup+self.outputInterval):self.outputInterval:finalTime).';
            end
            t(t<initialTime) = [];
            t(t<self.initialTime) = [];
            t(t>self.finalTime) = [];
        end

    end

    methods (Static)
        function vars = classRequiredPropertyNames()
            vars = {'observingSystems','name','outputInterval','initialTime','finalTime'};
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            arguments (Output)
                propertyAnnotations CAPropertyAnnotation
            end
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
            propertyAnnotations(end+1) = CAObjectProperty('observingSystems','array of WVObservingSystem objects');
            propertyAnnotations(end+1) = CAPropertyAnnotation('name','name of output group');
            propertyAnnotations(end+1) = CANumericProperty('outputInterval', {}, 's','output interval');
            propertyAnnotations(end+1) = CANumericProperty('initialTime', {}, 's','model time of first allowed write to file');
            propertyAnnotations(end+1) = CANumericProperty('finalTime', {}, 's','model time of last allowed write to file');
        end
    end
end