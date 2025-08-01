classdef WVModelOutputGroup < handle & matlab.mixin.Heterogeneous & CAAnnotatedClass
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
        function self = WVModelOutputGroup(model,options)
            arguments
                model WVModel
                options.name {mustBeText}
            end
            if ~isfield(options,"name")
                error("You must specify an output group name");
            end
            self.model = model;
            self.name = options.name;
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
                observingSystem WVObservingSystem
            end
            if self.didInitializeStorage
                error('Storage already initialized! You cannot add a new observing system after the storage has been initialized.');
            end

            for iObs = 1:length(observingSystem)
                anObservingSystem = observingSystem(iObs);
                if anObservingSystem.wvt ~= self.model.wvt
                    error('This observing system was not initialized with the same wvt that it is being added to!')
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
                observingSystem WVObservingSystem
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

        function initializeOutputGroup(self,ncfile)
            arguments (Input)
                self WVModelOutputGroup {mustBeNonempty}
                ncfile NetCDFGroup {mustBeNonempty}
            end
            if self.didInitializeStorage
                error('Storage already initialized!');
            end
            self.group = ncfile.addGroup(self.name);
            self.writeToGroup(self.group,self.propertyAnnotationWithName(self.requiredProperties));

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
            self.didInitializeStorage = true;
        end

        function writeTimeStepToNetCDFFile(self,ncfile,t)
            arguments
                self WVModelOutputGroup
                ncfile NetCDFFile
                t double
            end
            if ~self.didInitializeStorage
                self.initializeOutputGroup(ncfile);
            end
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

        function recordNetCDFFileHistory(self,options)
            arguments
                self WVModelOutputGroup {mustBeNonempty}
                options.didBlowUp {mustBeNumeric} = 0
            end
            if isempty(self.group)
                return
            end

            if options.didBlowUp == 1
                a = sprintf('%s: wrote %d time points to file. Terminated due to model blow-up.',datetime('now'),self.incrementsWrittenToGroup);
            else
                a = sprintf('%s: wrote %d time points to file',datetime('now'),self.incrementsWrittenToGroup);
            end
            if isKey(self.group.attributes,'history')
                history = reshape(self.group.attributes('history'),1,[]);
                history =cat(2,squeeze(history),a);
            else
                history = a;
            end
            self.group.addAttribute('history',history);
        end

        function closeNetCDFFile(self)
            if ~isempty(self.group)
                fprintf('Ending simulation. Wrote %d time points to %s group\n',self.incrementsWrittenToGroup,self.name);
            end
        end

        function initObservingSystemsFromGroup(self,outputGroup)
            arguments
                self WVModelOutputGroup {mustBeNonempty}
                outputGroup NetCDFGroup {mustBeNonempty}
            end

            f = @(className,group) feval(strcat(className,'.observingSystemFromGroup'),group, self.model, self);
            vars = CAAnnotatedClass.propertyValuesFromGroup(outputGroup,{"observingSystems"},classConstructor=f);
            self.addObservingSystem(vars.observingSystems);
        end

    end

    methods (Static)
        function outputGroup = modelOutputGroupFromGroup(group,model)
            %initialize a WVModelOutputGroup instance from NetCDF file
            %
            % Subclasses to should override this method to enable model
            % restarts. This method works in conjunction with -writeToFile
            % to provide restart capability.
            %
            % - Topic: Initialization
            % - Declaration: outputGroup = modelOutputGroupFromGroup(group,wvt)
            % - Parameter wvt: the WVTransform to be used
            % - Returns outputGroup: a new instance of WVModelOutputGroup
            arguments
                group NetCDFGroup {mustBeNonempty}
                model WVModel {mustBeNonempty}
            end
            className = group.attributes('AnnotatedClass');
            requiredProperties = feval(strcat(className,'.classRequiredPropertyNames'));
            requiredProperties(ismember(requiredProperties,'observingSystems')) = [];
            vars = CAAnnotatedClass.propertyValuesFromGroup(group,requiredProperties);
            if isempty(vars)
                outputGroup = feval(className,model);
            else
                options = namedargs2cell(vars);
                outputGroup = feval(className,model,options{:});
            end
            outputGroup.group = group;
            nPoints = group.dimensionWithName("t").nPoints;
            outputGroup.incrementsWrittenToGroup = nPoints;
            if nPoints > 0
                outputGroup.timeOfLastIncrementWrittenToGroup = group.readVariablesAtIndexAlongDimension('t',nPoints,'t');
            end
            outputGroup.initObservingSystemsFromGroup(group);
            outputGroup.didInitializeStorage = true;
        end
    end
end