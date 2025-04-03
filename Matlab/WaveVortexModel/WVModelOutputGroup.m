classdef WVModelOutputGroup < handle & matlab.mixin.Heterogeneous
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties (WeakHandle)
        % Reference to the WVModel being used
        model WVModel

        % Reference to the NetCDFGroup being used for model output
        % - Topic: Writing to NetCDF files
        % Empty indicates no file output. The output group creates the
        % NetCDFGroup, but the NetCDFFile owns it, hence a WeakHandle.
        group NetCDFGroup
    end

    properties
        name string

        % output index of the current/most recent step.
        % - Topic: Integration
        % If stepsTaken=0, outputIndex=1 means the initial conditions get written at index 1
        incrementsWrittenToGroup (1,1) uint64 = 0

        % output index of the current/most recent step.
        % - Topic: Integration
        % If stepsTaken=0, outputIndex=1 means the initial conditions get written at index 1
        timeOfLastIncrementWrittenToGroup (1,1) double = -Inf

        didInitializeStorage = false
    end

    properties
        observingSystems
    end

    methods
        function self = WVModelOutputGroup(model,name)
            arguments
                model WVModel
                name {mustBeText}
            end
            self.model = model;
            self.name = name;
            self.observingSystems = WVObservingSystem.empty(1,0);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Add/remove observing systems
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function addObservingSystem(self,observingSystem)
            arguments
                self WVModelOutputGroup {mustBeNonempty}
                observingSystem WVModelOutputGroup
            end
            if self.didInitializeStorage
                error('Storage already initialized! You cannot add a new observing system after the storage has been initialized.');
            end

            for iObs = 1:length(observingSystem)
                anObservingSystem = observingSystem(iObs);
                if anObservingSystem.wvt ~= self
                    error('This force was not initialized with the same wvt that it is being added to!')
                end

                if ~isempty(find(strcmp({self.observingSystems.name}, anObservingSystem.name), 1))
                    error('An observing system named %s already exists.\n',anObservingSystem.name)
                end

                self.observingSystems(end+1) = anObservingSystem;
                if anObservingSystem.nFluxComponents > 0
                    self.model.addFluxedObservingSystem(anObservingSystem);
                end
            end
        end

        function removeObservingSystem(self, observingSystem)
            arguments
                self WVModelOutputGroup {mustBeNonempty}
                observingSystem WVModelOutputGroup
            end

            for iObs = 1:length(observingSystem)
                anObservingSystem = observingSystem(iObs);

                % Verify that the observing system belongs to the same wvt
                if anObservingSystem.wvt ~= self
                    error('This observing system does not belong to the same wvt!');
                end

                % Find index of the observing system with the matching name
                idx = find(strcmp({self.observingSystems.name}, anObservingSystem.name));
                if isempty(idx)
                    error('No observing system named %s exists to remove.', anObservingSystem.name);
                end

                % Remove the observing system from the collection
                self.observingSystems(idx) = [];

                % If the system includes flux components, remove it from the fluxed systems in the model
                if anObservingSystem.nFluxComponents > 0
                    self.model.removeFluxedObservingSystem(anObservingSystem);
                end
            end
        end

        function observingSystem = observingSystemWithName(self,name)
            idx = find(strcmp({self.observingSystems.name}, name),1);
            if isempty(idx)
                error('No observing system named %s exists to remove.', anObservingSystem.name);
            else
                observingSystem = self.observingSystems(idx);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Write to file
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function t = outputTimesForIntegrationPeriod(self,initialTime,finalTime)
            arguments (Input)
                self WVModelOutputGroup
                initialTime (1,1) double
                finalTime (1,1) double
            end
            arguments (Output)
                t (:,1) double
            end
            t = [];
        end

        function openNetCDFFileForTimeStepping(self,ncfile)
            arguments (Input)
                self WVModelOutputGroup {mustBeNonempty}
                ncfile NetCDFGroup {mustBeNonempty}
            end
            if self.didInitializeStorage
                error('Storage already initialized!');
            end
            self.group = ncfile.addGroup(self.name);

            varAnnotation = self.model.wvt.propertyAnnotationWithName('t');
            varAnnotation.attributes('units') = varAnnotation.units;
            varAnnotation.attributes('long_name') = varAnnotation.description;
            varAnnotation.attributes('standard_name') = 'time';
            varAnnotation.attributes('long_name') = 'time';
            varAnnotation.attributes('units') = 'seconds since 1970-01-01 00:00:00';
            varAnnotation.attributes('axis') = 'T';
            varAnnotation.attributes('calendar') = 'standard';
            self.group.addDimension(varAnnotation.name,length=Inf,type="double",attributes=varAnnotation.attributes);

            for iObs = 1:length(self.observingSystems)
                self.observingSystems(iObs).initializeStorage(self.group);
            end

            self.incrementsWrittenToGroup = 0;
            self.writeTimeStepToNetCDFFile(self.model.t);
        end

        function writeTimeStepToNetCDFFile(self,t)
            if ( ~isempty(self.group) && t > self.timeOfLastIncrementWrittenToGroup )
                outputIndex = self.incrementsWrittenToGroup + 1;

                self.group.variableWithName('t').setValueAlongDimensionAtIndex(t,'t',outputIndex);

                for iObs = 1:length(self.observingSystems)
                    self.observingSystems(iObs).writeTimeStepToFile(self.group,outputIndex);
                end

                self.incrementsWrittenToGroup = outputIndex;
                self.timeOfLastIncrementWrittenToGroup = t;
            end
        end

        function closeNetCDFFile(self)
            if ~isempty(self.group)
                fprintf('Ending simulation. Wrote %d time points to %s group\n',self.incrementsWrittenToGroup,self.name);
            end
        end

    end
end