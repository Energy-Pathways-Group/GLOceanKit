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
        wvt
    end

    properties (Access=private)
        outputGroupNameMap = configureDictionary("string","WVModelOutputGroup")
        outputGroupNameOutputTimeMap = configureDictionary("string","cell")
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
            filename = strcat(name,ext);
        end

        function wvt = get.wvt(self)
            wvt = self.model.wvt;
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

        function addOutputGroup(self,outputGroup)
            arguments
                self WVModelOutputFile {mustBeNonempty}
                outputGroup WVModelOutputGroup
            end
            self.outputGroupNameMap(outputGroup.name) = outputGroup;
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
                self.outputGroupNameOutputTimeMap{outputGroups_(iGroup).name} = t_group;
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
            outputGroupNames = self.outputGroupNameOutputTimeMap.keys;
            didWriteToFile = false;
            for i = 1:length(outputGroupNames)
                t_group = self.outputGroupNameOutputTimeMap{outputGroupNames(i)};
                if ~isempty(t_group) && abs(t - t_group(1)) < eps
                    self.outputGroupWithName(outputGroupNames(i)).writeTimeStepToNetCDFFile(t);
                    t_group(1) = [];
                    self.outputGroupNameOutputTimeMap{outputGroupNames(i)} = t_group;
                    didWriteToFile = true;
                end
            end
            if didWriteToFile
                self.ncfile.sync();
            end
        end

        function initializeOutputFile(self)
            if self.observingSystemWillWriteWaveVortexCoefficients == true
                properties = setdiff(self.wvt.requiredProperties,{'Ap','Am','A0','t'});
            else
                properties = setdiff(self.wvt.requiredProperties,{'t'});
            end
            % in theory we already removed the file if the user requested
            if isfile(self.path)
                error('A file already exists at this path.');
            end
            self.ncfile = self.wvt.writeToFile(self.path,properties{:},shouldOverwriteExisting=false,shouldAddRequiredProperties=false);
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
                        bool = bool | all(ismember(intersect({'Ap','Am','A0'},self.wvt.variableNames),observingSystems(iObs).netCDFOutputVariables));
                    end
                end
            end
        end

        function recordNetCDFFileHistory(self,options)
            arguments
                self WVModelOutputFile {mustBeNonempty}
                options.didBlowUp {mustBeNumeric} = 0
            end

            if ~isempty(self.ncfile)
                arrayfun( @(outputGroup) outputGroup.recordNetCDFFileHistory(didBlowUp=options.didBlowUp), self.outputGroups);
                self.ncfile = [];
            end
        end


        function closeNetCDFFile(self)
            if ~isempty(self.ncfile)
                arrayfun( @(outputGroup) outputGroup.closeNetCDFFile(), self.outputGroups);
                self.ncfile = [];
            end
        end

    end
end