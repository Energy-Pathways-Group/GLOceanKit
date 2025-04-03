classdef WVModelOutputFile < handle & matlab.mixin.Heterogeneous
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties (WeakHandle)
        % Reference to the WVModel being used
        model WVModel
    end

    properties
        % Reference to the NetCDFFile being used for model output
        ncfile NetCDFFile

        didInitializeStorage = false
        tInitialize = Inf;
    end

    properties (Dependent)
        outputGroups
    end

    properties (Access=private)
        outputGroupNameMap = configureDictionary("string","WVModelOutputGroup")
        outputGroupOutputTimeMap = configureDictionary("WVModelOutputGroup","double")
    end

    methods
        function self = WVModelOutputFile(model,name)
            arguments
                model WVModel
                name {mustBeText}
            end
            self.model = model;
            self.name = name;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Add/remove output groups
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function outputGroups = get.outputGroups(self)
            outputGroups = [self.outputGroupNameMap(self.outputGroupNameMap.keys)];
        end

        function names = outputGroupNames(self)
            % retrieve the names of all output group names
            %
            % - Topic: Utility function — Metadata
            arguments (Input)
                self WVModelOutputFile {mustBeNonempty}
            end
            arguments (Output)
                names string
            end
            names = self.outputGroupNameMap.keys;
        end

        function val = outputGroupWithName(self,name)
            % retrieve a WVModelOutputGroup by name
            arguments (Input)
                self WVModelOutputFile {mustBeNonempty}
                name char {mustBeNonempty}
            end
            arguments (Output)
                val WVModelOutputGroup
            end
            val = self.outputGroupNameMap(name);
        end

        function outputGroup = addOutputGroup(self,name,options)
            arguments
                self WVModelOutputFile {mustBeNonempty}
                name {mustBeText}
                options.outputInterval
            end
            if isfield(options,"outputInterval")
                outputGroup = WVModelOutputGroup(self,name,outputInterval=options.outputInterval);
            else
                outputGroup = WVModelOutputGroup(self,name);
            end
            self.outputGroupNameMap(name) = outputGroup;
        end

        function outputGroup = defaultOutputGroup(self)
            arguments (Input)
                self WVModelOutputFile {mustBeNonempty}
            end
            arguments (Output)
                outputGroup WVModelOutputGroup
            end
            outputGroup = self.outputGroupWithName(self.defaultOutputGroupName);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Write to file
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function t = outputTimesForIntegrationPeriod(self,initialTime,finalTime)
            % This will be called exactly once before an integration
            % begins.
            arguments (Input)
                self WVModelOutputFile
                initialTime (1,1) double
                finalTime (1,1) double
            end
            arguments (Output)
                t (:,1) double
            end
            t = [];
            outputGroups_ = self.outputGroups;
            for iGroup = 1:length(outputGroups_)
                t_group = outputGroups_(iGroup).outputTimesForIntegrationPeriod(initialTime,finalTime);
                self.outputGroupOutputTimeMap(outputGroups_(iGroup)) = t_group;
                t = cat(1,t,t_group);
            end
            t = sort(t);
            if self.didInitializeStorage == false && ~isempty(t)
                self.tInitialize = t(1);
            end
        end

        function writeTimeStepToOutputFile(self,t)
            % Here we need to insert the logical to create the NetCDF file
            % if it hasn't been created.
            % 
            outputGroups_ = self.outputGroupOutputTimeMap.keys;
            for i = 1:length(outputGroups_)
                t_group = self.outputGroupOutputTimeMap(outputGroups_(i));
                if ~isempty(t_group) && abs(t - t_group(1)) < eps
                    outputGroups_(i).writeTimeStepToNetCDFFile(t);
                    t_group(1) = [];
                    self.outputGroupOutputTimeMap(outputGroups_(i)) = t_group;
                end
            end
        end

        function bool = observingSystemWillWriteWaveVortexCoefficients(self)
            % A simple check to see if one of the observing systems will be
            % writing wave-vortex coefficients
            outputGroups_ = self.outputGroups;
            bool = false;
            for iGroup = 1:length(outputGroups_)
                observingSystems = outputGroups_(iGroup).observingSystems;
                for iObs = 1:length(observingSystems)
                    if isa(observingSystems(iObs),'WVEulerianFields')
                        bool = bool | all(ismember(intersect({'Ap','Am','A0'},wvt.variableNames),observingSystems(iObs).netCDFOutputVariables));
                    end
                end
            end
        end

        function openNetCDFFileForTimeStepping(self,ncfile)
            arguments (Input)
                self WVModelOutputFile {mustBeNonempty}
                ncfile NetCDFGroup {mustBeNonempty}
            end
            if self.didInitializeStorage
                error('Storage already initialized!');
            end

% abs(t(iTime) - self.outputTimes(1)) < eps

            self.ncfile = ncfile.addGroup(self.name);

            varAnnotation = self.model.wvt.propertyAnnotationWithName('t');
            varAnnotation.attributes('units') = varAnnotation.units;
            varAnnotation.attributes('long_name') = varAnnotation.description;
            varAnnotation.attributes('standard_name') = 'time';
            varAnnotation.attributes('long_name') = 'time';
            varAnnotation.attributes('units') = 'seconds since 1970-01-01 00:00:00';
            varAnnotation.attributes('axis') = 'T';
            varAnnotation.attributes('calendar') = 'standard';
            self.ncfile.addDimension(varAnnotation.name,length=Inf,type="double",attributes=varAnnotation.attributes);

            for iObs = 1:length(self.observingSystems)
                self.observingSystems(iObs).initializeStorage(self.ncfile);
            end

            self.incrementsWrittenToGroup = 0;
            self.writeTimeStepToNetCDFFile(self.model.t);
        end

        function writeTimeStepToNetCDFFile(self,t)
            if ( ~isempty(self.ncfile) && t > self.timeOfLastIncrementWrittenToGroup )
                outputIndex = self.incrementsWrittenToGroup + 1;

                self.ncfile.variableWithName('t').setValueAlongDimensionAtIndex(t,'t',outputIndex);

                for iObs = 1:length(self.observingSystems)
                    self.observingSystems(iObs).writeTimeStepToFile(self.ncfile,outputIndex);
                end

                self.incrementsWrittenToGroup = outputIndex;
                self.timeOfLastIncrementWrittenToGroup = t;
            end
        end

        function closeNetCDFFile(self)
            if ~isempty(self.ncfile)
                fprintf('Ending simulation. Wrote %d time points to %s group\n',self.incrementsWrittenToGroup,self.name);
            end
        end

    end
end