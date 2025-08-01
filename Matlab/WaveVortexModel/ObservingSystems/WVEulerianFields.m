classdef WVEulerianFields < WVObservingSystem
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties (GetAccess=public, SetAccess=protected)
        initialConditionOnlyVariables = {};
        timeSeriesVariables = {};
    end

    properties (SetObservable)
        netCDFOutputVariables = {};
    end

    properties (Dependent)
        nOutputVariables
        nTimeSeriesVariables
        fieldNames
    end

    methods
        function self = WVEulerianFields(model,options)
            %create a new observing system
            %
            % This class is intended to be subclassed, so it generally
            % assumed that this initialization will not be called directly.
            %
            % - Topic: Initialization
            % - Declaration: self = WVObservingSystem(model,name)
            % - Parameter model: the WVModel instance
            % - Parameter name: name of the observing system
            % - Returns self: a new instance of WVObservingSystem
            arguments
                model WVModel
                options.fieldNames
            end
            % Do we actually want to inherit the properties from the
            % WVTransform? I'm not sure. I think this should be optional.
            % If an OS does, then its output can go in the wave-vortex
            % group.
            if ~isfield(options,"fieldNames")
                options.fieldNames = {};
            elseif isa(options.fieldNames,"string")
                options.fieldNames = cellstr(options.fieldNames);
            end
            self@WVObservingSystem(model,"eulerian fields");
            addlistener(self,'netCDFOutputVariables','PostSet',@(src,evnt) self.updateNetCDFVariableCategorization);
            self.netCDFOutputVariables = options.fieldNames;
        end

        function names = get.fieldNames(self)
            names = string(self.netCDFOutputVariables);
        end

        function nOutputVariables = get.nOutputVariables(self)
            nOutputVariables = length(self.netCDFOutputVariables);
        end

        function nOutputVariables = get.nTimeSeriesVariables(self)
            nOutputVariables = 0;
            for iVar = 1:length(self.netCDFOutputVariables)
                varAnnotation = self.model.wvt.propertyAnnotationWithName(self.netCDFOutputVariables{iVar});
                if (self.model.isDynamicsLinear == 1 && varAnnotation.isVariableWithLinearTimeStep == 1) || (self.model.isDynamicsLinear == 0 && varAnnotation.isVariableWithNonlinearTimeStep == 1)
                    nOutputVariables = nOutputVariables + 1;
                end
            end
        end

        function aString = description(self)
            aString = description@WVObservingSystem(self);
            aString = aString + ", writing fields " + strjoin(string(self.timeSeriesVariables),", ");
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Add/remove output variables
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
                self WVEulerianFields
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
                self WVEulerianFields
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
                self WVEulerianFields
            end
            arguments (Repeating)
                variables char
            end
            self.netCDFOutputVariables = setdiff(self.netCDFOutputVariables,variables);
        end

        function updateNetCDFVariableCategorization(self)
            self.initialConditionOnlyVariables = {};
            self.timeSeriesVariables = {};
            for iVar = 1:length(self.netCDFOutputVariables)
                varAnnotation = self.model.wvt.propertyAnnotationWithName(self.netCDFOutputVariables{iVar});
                if (self.model.isDynamicsLinear == 1 && varAnnotation.isVariableWithLinearTimeStep == 1) || (self.model.isDynamicsLinear == 0 && varAnnotation.isVariableWithNonlinearTimeStep == 1)
                    self.timeSeriesVariables{end+1} = self.netCDFOutputVariables{iVar};
                else
                    self.initialConditionOnlyVariables{end+1} = self.netCDFOutputVariables{iVar};
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Read and write to file
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function initializeStorage(self,group)
            % Sort through which variables we will record a time series
            % for, and which we will only write initial conditions.
            self.initialConditionOnlyVariables = {};
            self.timeSeriesVariables = {};
            for iVar = 1:length(self.netCDFOutputVariables)
                if group.hasVariableWithName(self.netCDFOutputVariables{iVar})
                    continue;
                end

                varAnnotation = self.model.wvt.propertyAnnotationWithName(self.netCDFOutputVariables{iVar});
                varAnnotation.attributes('units') = varAnnotation.units;
                varAnnotation.attributes('long_name') = varAnnotation.description;

                for iDim=1:length(varAnnotation.dimensions)
                    if ~group.hasDimensionWithName(varAnnotation.dimensions{iDim})
                        dimAnnotation = self.model.wvt.propertyAnnotationWithName(varAnnotation.dimensions{iDim});
                        dimAnnotation.attributes('units') = dimAnnotation.units;
                        dimAnnotation.attributes('long_name') = dimAnnotation.description;
                        group.addDimension(dimAnnotation.name,self.model.wvt.(dimAnnotation.name),attributes=dimAnnotation.attributes);
                    end
                end

                if (self.model.isDynamicsLinear == 1 && varAnnotation.isVariableWithLinearTimeStep == 1) || (self.model.isDynamicsLinear == 0 && varAnnotation.isVariableWithNonlinearTimeStep == 1)
                    self.timeSeriesVariables{end+1} = self.netCDFOutputVariables{iVar};
                    group.addVariable(varAnnotation.name,horzcat(varAnnotation.dimensions,'t'),type="double",isComplex=varAnnotation.isComplex,attributes=varAnnotation.attributes);
                else
                    self.initialConditionOnlyVariables{end+1} = self.netCDFOutputVariables{iVar};
                    group.addVariable(varAnnotation.name,varAnnotation.dimensions,self.model.wvt.(varAnnotation.name),attributes=varAnnotation.attributes);
                end
            end
        end

        function writeTimeStepToFile(self,group,outputIndex)
            for iVar=1:length(self.timeSeriesVariables)
                group.variableWithName(self.timeSeriesVariables{iVar}).setValueAlongDimensionAtIndex(self.model.wvt.(self.timeSeriesVariables{iVar}),'t',outputIndex);
            end
        end

        function os = observingSystemWithResolutionOfTransform(self,wvtX2)
            %create a new WVObservingSystem with a new resolution
            %
            % Subclasses to should override this method an implement the
            % correct logic.
            %
            % - Topic: Initialization
            % - Declaration: os = observingSystemWithResolutionOfTransform(self,wvtX2)
            % - Parameter wvtX2: the WVTransform with increased resolution
            % - Returns force: a new instance of WVObservingSystem
            os = WVEulerianFields(wvtX2,self.name);
        end
    end

    methods (Static)
        function vars = classRequiredPropertyNames()
            vars = {'fieldNames'};
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            arguments (Output)
                propertyAnnotations CAPropertyAnnotation
            end
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
            propertyAnnotations(end+1) = CAPropertyAnnotation('fieldNames','eulerian field names');
        end
    end
end