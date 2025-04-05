classdef WVModelOutputFile < handle & matlab.mixin.Heterogeneous
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties (WeakHandle)
        % Reference to the WVModel being used
        model WVModel
    end

    properties
        path

        % Reference to the NetCDFFile being used for model output
        ncfile NetCDFFile

        didInitializeStorage = false
        tInitialize = Inf;
    end

    properties (Dependent)
        outputGroups
        filename
    end

    properties (Access=private)
        outputGroupNameMap = configureDictionary("string","WVModelOutputGroup")
        outputGroupOutputTimeMap = configureDictionary("WVModelOutputGroup","double")
    end

    methods
        function self = WVModelOutputFile(model,path,options)
            arguments
                model WVModel
                path {mustBeText}
                options.shouldOverwriteExisting double {mustBeMember(options.shouldOverwriteExisting,[0 1])} = 0
            end
            if options.shouldOverwriteExisting == 1
                if isfile(path)
                    delete(path);
                end
            else
                if isfile(path)
                    error('A file already exists with that name.')
                end
            end
            self.model = model;
            self.path = path;
        end

        function filename = get.filename(self)
            [~,name,ext] = fileparts(self.path);
            filename = name + ext;
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
            % 1) initialize the netcdf file if necessary
            if self.didInitializeStorage == false && abs(t - self.tInitialize) < eps
                self.initializeOutputFile();
            end

            % 2) inform the appropriate groups that they need to write a time step.
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

        function initializeOutputFile(self)
            if self.observingSystemWillWriteWaveVortexCoefficients == true
                properties = setdiff(self.wvt.requiredProperties,{'Ap','Am','A0','t'});
            else
                properties = setdiff(self.wvt.requiredProperties,{'t'});
            end
            self.ncfile = self.wvt.writeToFile(netcdfFile,properties{:},shouldOverwriteExisting=options.shouldOverwriteExisting,shouldAddRequiredProperties=false);
            self.didInitializeStorage = true;

            arrayfun( @(outputGroup) outputGroup.initializeOutputGroup(self.ncfile), self.outputGroups);
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

        function recordNetCDFFileHistory(self,options)
            arguments
                self WVModelOutputFile {mustBeNonempty}
                options.didBlowUp {mustBeNumeric} = 0
            end
            if isempty(self.ncfile)
                return
            end

            if options.didBlowUp == 1
                a = sprintf('%s: wrote %d time points to file. Terminated due to model blow-up.',datetime('now'),self.incrementsWrittenToFile);
            else
                a = sprintf('%s: wrote %d time points to file',datetime('now'),self.incrementsWrittenToFile);
            end
            if isKey(self.ncfile.attributes,'history')
                history = reshape(self.ncfile.attributes('history'),1,[]);
                history =cat(2,squeeze(history),a);
            else
                history = a;
            end
            self.ncfile.addAttribute('history',history);
        end


        function closeNetCDFFile(self)
            if ~isempty(self.ncfile)
                arrayfun( @(outputGroup) outputGroup.closeNetCDFFile(), self.outputGroups);
                self.ncfile = [];
            end
        end

    end
end