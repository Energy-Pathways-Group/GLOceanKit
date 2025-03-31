classdef WVModelOutputGroup < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        model WVModel

        name string

        % Reference to the NetCDFFile being used for model output
        % - Topic: Writing to NetCDF files
        % Empty indicates no file output.
        group NetCDFGroup

        % List of all variables being written to NetCDF file
        % - Topic: Writing to NetCDF files
        % The default list includes 'Ap', 'Am', and 'A0' which is the
        % minimum set of variables required for a restart.
        netCDFOutputVariables = {}

        netCDFOutputParticles = {};

        netCDFOutputTracers = {};

        % Model output interval (seconds)
        % - Topic: Integration
        % This property is optionally set when calling setupIntegrator. If
        % set, it will allow you to call -integrateToNextOutputTime and, if
        % a NetCDF file is set for output, it will set the interval at
        % which time steps are written to file.
        outputInterval = [] % (1,1) double

        % output index of the current/most recent step.
        % - Topic: Integration
        % If stepsTaken=0, outputIndex=1 means the initial conditions get written at index 1
        incrementsWrittenToFile (1,1) uint64 = 0

        % output index of the current/most recent step.
        % - Topic: Integration
        % If stepsTaken=0, outputIndex=1 means the initial conditions get written at index 1
        timeOfLastIncrementWrittenToFile (1,1) double = -Inf
    end

    properties
        % map to a map containing the particle variables, e.g. particlesWithName('float') returns a map containing keys ('x','y','z') at minimum
        netcdfVariableMapForParticleWithName
        initialConditionOnlyVariables
        timeSeriesVariables
    end

    methods
        function self = WVModelOutputGroup(model,name,options)
            arguments
                model WVModel
                name {mustBeText}
                options.outputInterval (1,1) double {mustBePositive}
            end
            self.model = model;
            self.name = name;
            self.netcdfVariableMapForParticleWithName = containers.Map();
            if isfield(options,"outputInterval")
                self.outputInterval = options.outputInterval;
            end
        end

        function addNetCDFOutputParticles(self,particles)
            self.netCDFOutputParticles = union(self.netCDFOutputParticles,particles);
        end
        
        function addNetCDFOutputTracers(self,tracers)
            self.netCDFOutputTracers = union(self.netCDFOutputTracers,tracers);
        end

        function addNetCDFOutputVariables(self,variables)
            % Add variables to list of variables to be written to the NetCDF variable during the model run.
            %
            % - Topic: Writing to NetCDF files
            % - Declaration: addNetCDFOutputVariables(variables)
            % - Parameter variables: strings of variable names.
            %
            % Pass strings of WVTransform state variables of the
            % same name. This must be called before using any of the
            % integrate methods.
            %
            % ```matlab
            % model.addNetCDFOutputVariables('A0','u','v');
            % ```

            arguments
                self WVModelOutputGroup
            end
            arguments (Repeating)
                variables char
            end
            unknownVars = setdiff(variables,self.model.wvt.variableNames);
            if ~isempty(unknownVars)
                error('The WVTransform does not have a variable named %s',unknownVars{1}) ;
            end
            self.netCDFOutputVariables = union(self.netCDFOutputVariables,variables);
        end

        function setNetCDFOutputVariables(self,variables)
            % Set list of variables to be written to the NetCDF variable during the model run.
            %
            % - Topic: Writing to NetCDF files
            % - Declaration: setNetCDFOutputVariables(variables)
            % - Parameter variables: strings of variable names.
            %
            % Pass strings of WVTransform state variables of the
            % same name. This must be called before using any of the
            % integrate methods.
            %
            % ```matlab
            % model.setNetCDFOutputVariables('A0','u','v');
            % ```
            arguments
                self WVModelOutputGroup
            end
            arguments (Repeating)
                variables char
            end
            unknownVars = setdiff(variables,self.model.wvt.variableNames);
            if ~isempty(unknownVars)
                error('The WVTransform does not have a variable named %s',unknownVars{1}) ;
            end
            self.netCDFOutputVariables = variables;
        end

        function removeNetCDFOutputVariables(self,variables)
            % Remove variables from the list of variables to be written to the NetCDF variable during the model run.
            %
            % - Topic: Writing to NetCDF files
            % - Declaration: removeNetCDFOutputVariables(variables)
            % - Parameter variables: strings of variable names.
            %
            % Pass strings of WVTransform state variables of the
            % same name. This must be called before using any of the
            % integrate methods.
            %
            % ```matlab
            % model.removeNetCDFOutputVariables('A0','u','v');
            % ```
            arguments
                self WVModelOutputGroup
            end
            arguments (Repeating)
                variables char
            end
            self.netCDFOutputVariables = setdiff(self.netCDFOutputVariables,variables);
        end

        function openNetCDFFileForTimeStepping(self,ncfile)
            arguments (Input)
                self WVModelOutputGroup {mustBeNonempty}
                ncfile NetCDFGroup {mustBeNonempty}
            end
            % if isempty(self.group)
                if isempty(self.outputInterval)
                    error('The output interval for this group was not set!!!')
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

                % Sort through which variables we will record a time series
                % for, and which we will only write initial conditions.
                self.initialConditionOnlyVariables = {};
                self.timeSeriesVariables = {};
                for iVar = 1:length(self.netCDFOutputVariables)
                    if self.group.hasVariableWithName(self.netCDFOutputVariables{iVar})
                        continue;
                    end

                    varAnnotation = self.model.wvt.propertyAnnotationWithName(self.netCDFOutputVariables{iVar});
                    varAnnotation.attributes('units') = varAnnotation.units;
                    varAnnotation.attributes('long_name') = varAnnotation.description;

                    for iDim=1:length(varAnnotation.dimensions)
                        if ~self.group.hasDimensionWithName(varAnnotation.dimensions{iDim})
                            dimAnnotation = self.model.wvt.propertyAnnotationWithName(varAnnotation.dimensions{iDim});
                            dimAnnotation.attributes('units') = dimAnnotation.units;
                            dimAnnotation.attributes('long_name') = dimAnnotation.description;
                            self.group.addDimension(dimAnnotation.name,self.model.wvt.(dimAnnotation.name),attributes=dimAnnotation.attributes);
                        end
                    end

                    if (self.model.linearDynamics == 1 && varAnnotation.isVariableWithLinearTimeStep == 1) || (self.model.linearDynamics == 0 && varAnnotation.isVariableWithNonlinearTimeStep == 1)
                        self.timeSeriesVariables{end+1} = self.netCDFOutputVariables{iVar};
                        self.group.addVariable(varAnnotation.name,horzcat(varAnnotation.dimensions,'t'),type="double",isComplex=varAnnotation.isComplex,attributes=varAnnotation.attributes);
                    else
                        self.initialConditionOnlyVariables{end+1} = self.netCDFOutputVariables{iVar};
                        self.group.addVariable(varAnnotation.name,varAnnotation.dimensions,self.model.wvt.(varAnnotation.name),attributes=varAnnotation.attributes);
                    end
                end

                for iTracer = 1:length(self.netCDFOutputTracers)
                    if self.group.hasVariableWithName(self.netCDFOutputTracers{iTracer})
                        continue;
                    end
                    self.group.addVariable(self.netCDFOutputTracers{iTracer},horzcat(self.model.wvt.spatialDimensionNames,'t'),type="double",attributes=containers.Map({'isTracer'},{true}));
                end

                for iParticle = 1:length(self.netCDFOutputParticles)
                    if self.group.hasVariableWithName(self.netCDFOutputParticles{iParticle})
                        continue;
                    end
                    particle = self.model.particle{self.model.particleIndexWithName(self.netCDFOutputParticles{iParticle})};
                    self.initializeParticleStorage(particle.name,size(particle.x,2),particle.trackedFieldNames{:});
                end
                self.incrementsWrittenToFile = 0;
                self.writeTimeStepToNetCDFFile(self.model.t);
            % end

        end

        function writeTimeStepToNetCDFFile(self,t)
            if ( ~isempty(self.group) && t > self.timeOfLastIncrementWrittenToFile )
                outputIndex = self.incrementsWrittenToFile + 1;

                self.group.variableWithName('t').setValueAlongDimensionAtIndex(t,'t',outputIndex);

                for iVar=1:length(self.timeSeriesVariables)
                    self.group.variableWithName(self.timeSeriesVariables{iVar}).setValueAlongDimensionAtIndex(self.model.wvt.(self.timeSeriesVariables{iVar}),'t',outputIndex);
                end

                for iParticle = 1:length(self.netCDFOutputParticles)
                    name_ = self.netCDFOutputParticles{iParticle};
                    [x,y,z,trackedFields] = self.model.particlePositions(name_);
                    self.writeParticleDataAtTimeIndex(name_,outputIndex,x,y,z,trackedFields);
                end

                for iTracer = 1:length(self.netCDFOutputTracers)
                    name_ = self.netCDFOutputTracers{iTracer};
                    self.group.variableWithName(name_).setValueAlongDimensionAtIndex(self.model.tracer(name_),'t',outputIndex);
                end

                self.incrementsWrittenToFile = outputIndex;
                self.timeOfLastIncrementWrittenToFile = t;
            end
        end

        function closeNetCDFFile(self)
            if ~isempty(self.group)
                fprintf('Ending simulation. Wrote %d time points to %s group\n',self.incrementsWrittenToFile,self.name);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Particles
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function initializeParticleStorage(self,particleName, nParticles, trackedFieldNames)
            arguments
                self WVModelOutputGroup
                particleName char
                nParticles (1,1) double {mustBePositive}
            end
            arguments (Repeating)
                trackedFieldNames char
            end

            variables = containers.Map();

            commonKeys = {'isParticle','particleName'};
            commonVals = {1,particleName};
            attributes = containers.Map(commonKeys,commonVals);
            attributes('units') = 'unitless id number';
            attributes('particleVariableName') = 'id';
            [dim,var] = self.group.addDimension(strcat(particleName,'_id'),(1:nParticles).',attributes=attributes);
            variables('id') = var;
            
            % careful to create a new object each time we init
            spatialDimensionNames = self.model.wvt.spatialDimensionNames;
            for iVar=1:length(spatialDimensionNames)
                attributes = containers.Map(commonKeys,commonVals);
                attributes('units') = self.model.wvt.dimensionAnnotationWithName(spatialDimensionNames{iVar}).units;
                attributes('long_name') = strcat(self.model.wvt.dimensionAnnotationWithName(spatialDimensionNames{iVar}).description,', recorded along the particle trajectory');
                attributes('particleVariableName') = spatialDimensionNames{iVar};
                variables(spatialDimensionNames{iVar}) = self.group.addVariable(strcat(particleName,'_',spatialDimensionNames{iVar}),{dim.name,'t'},type="double",attributes=attributes,isComplex=false);
            end

            for iVar=1:length(trackedFieldNames)
                varAnnotation = self.model.wvt.propertyAnnotationWithName(trackedFieldNames{iVar});
                attributes = containers.Map(commonKeys,commonVals);
                attributes('units') = varAnnotation.units;
                attributes('long_name') = strcat(varAnnotation.description,', recorded along the particle trajectory');
                attributes('particleVariableName') = trackedFieldNames{iVar};
                variables(trackedFieldNames{iVar}) = self.group.addVariable(strcat(particleName,'_',trackedFieldNames{iVar}),{dim.name,'t'},type="double",attributes=attributes,isComplex=false);
            end
 
            self.netcdfVariableMapForParticleWithName(particleName) = variables;
        end

        function writeParticleDataAtTimeIndex(self,particleName,iTime,x,y,z,trackedFields)
            self.group.variableWithName(strcat(particleName,'_x')).setValueAlongDimensionAtIndex(x,'t',iTime);
            self.group.variableWithName(strcat(particleName,'_y')).setValueAlongDimensionAtIndex(y,'t',iTime);
            if ~isequal(self.model.wvt.spatialDimensionNames,{'x','y'})
                self.group.variableWithName(strcat(particleName,'_z')).setValueAlongDimensionAtIndex(z,'t',iTime);
            end

            if ~isempty(trackedFields)
                trackedFieldNames = fieldnames(trackedFields);
                for iField=1:length(trackedFieldNames)
                    self.group.variableWithName(strcat(particleName,'_',trackedFieldNames{iField})).setValueAlongDimensionAtIndex(trackedFields.(trackedFieldNames{iField}),'t',iTime);
                end
            end
        end
    end
end