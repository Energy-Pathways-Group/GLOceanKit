classdef WVLagrangianParticles < WVObservingSystem
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties (GetAccess=public, SetAccess=protected)
        x, y, z
        isXYOnly
        advectionInterpolation
        trackedFieldNamesCell
        trackedFields
        trackedVarInterpolation
        absToleranceXY
        absToleranceZ
    end

    properties (Dependent)
        nParticles
        trackedFieldNames
    end

    methods
        function self = WVLagrangianParticles(model,options)
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
                options.name {mustBeText}
                options.x (1,:) double
                options.y (1,:) double
                options.z (1,:) double
                options.isXYOnly (1,1) logical = false
                options.trackedFieldNames
                options.advectionInterpolation char {mustBeMember(options.advectionInterpolation,["linear","spline","exact","finufft"])} = "linear"
                options.trackedVarInterpolation char {mustBeMember(options.trackedVarInterpolation,["linear","spline","exact","finufft"])} = "linear"
                options.absToleranceXY = 1e-1; % 100 km * 10^{-6}
                options.absToleranceZ = 1e-2;  
            end
            % Do we actually want to inherit the properties from the
            % WVTransform? I'm not sure. I think this should be optional.
            % If an OS does, then its output can go in the wave-vortex
            % group.
            self@WVObservingSystem(model,options.name);
            if ~isfield(options,"trackedFieldNames")
                options.trackedFieldNames = {};
            elseif isa(options.trackedFieldNames,"string")
                options.trackedFieldNames = cellstr(options.trackedFieldNames);
            end

            % Confirm that we really can track these variables.
            for iVar=1:length(options.trackedFieldNames)
                if ~any(ismember(self.wvt.variableNames,options.trackedFieldNames{iVar}))
                    error('Unable to find a WVVariableAnnotation named %s.', options.trackedFieldNames{iVar});
                end
                transformVar = self.wvt.propertyAnnotationWithName(options.trackedFieldNames{iVar});
                if isequal(self.wvt.spatialDimensionNames,{'x','y'})
                    if ~all(ismember(transformVar.dimensions,{'x','y'})) && ~all(ismember(transformVar.dimensions,{'x','y','z'}))
                        error('The WVVariableAnnotation %s does not have dimensions (x,y) or (x,y,z) and theforefore cannot be used for particle tracking', options.trackedFieldNames{iVar});
                    end
                else
                    if ~all(ismember(transformVar.dimensions,{'x','y','z'}))
                        error('The WVVariableAnnotation %s does not have dimensions x,y,z and theforefore cannot be used for particle tracking', options.trackedFieldNames{iVar});
                    end
                end
            end

            self.x = options.x;
            self.y = options.y;
            self.z = options.z;
            self.isXYOnly = options.isXYOnly;
            self.trackedFieldNamesCell = options.trackedFieldNames;
            self.trackedVarInterpolation = options.trackedVarInterpolation;
            self.advectionInterpolation = options.advectionInterpolation;
            self.absToleranceXY = options.absToleranceXY;
            self.absToleranceZ = options.absToleranceZ;

            if self.isXYOnly
                self.nFluxComponents = 2;
            else
                self.nFluxComponents = 3;
            end

            self.updateParticleTrackedFields();
        end

        function names = get.trackedFieldNames(self)
            names = string(self.trackedFieldNamesCell);
        end

        function nParticles = get.nParticles(self)
            nParticles = length(self.x);
        end

        function [x,y,z,trackedFields] = particlePositions(self)
            % Positions and values of tracked fields of particles at the current model time.
            %
            % - Topic: Particles
            % - Declaration: [x,y,z,trackedFields] = particlePositions(name)
            % - Parameter name: name of the particles
            x = self.x;
            y = self.y;
            z = self.z;
            trackedFields = self.trackedFields;
        end

        function updateParticleTrackedFields(self)
            % One special thing we have to do is log the particle
            % tracked fields
            if ~isempty(self.trackedFieldNamesCell)
                varLagrangianValues = cell(1,length(self.trackedFieldNamesCell));
                [varLagrangianValues{:}] = self.wvt.variableAtPositionWithName(self.x,self.y,self.z,self.trackedFieldNamesCell{:},interpolationMethod=self.trackedVarInterpolation);
                for i=1:length(self.trackedFieldNamesCell)
                    self.trackedFields.(self.trackedFieldNamesCell{i}) = varLagrangianValues{i};
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Integrated variables
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function nArray = lengthOfFluxComponents(self)
            % return an array containing the numel of each flux component.
            if self.isXYOnly
                nArray = [numel(self.x);numel(self.y)];
            else
                nArray = [numel(self.x);numel(self.y);numel(self.z)];
            end
        end

        function Y0 = absErrorTolerance(self)
            % return a cell array of the absolute tolerances of the
            % variables being integrated. You can pass either scalar
            % values, or an array of the same size as the variable.
            %
            % this will only be called when the time-stepping is run with
            % an adaptive integrator.
            if self.isXYOnly
                Y0 = {self.absToleranceXY,self.absToleranceXY};
            else
                Y0 = {self.absToleranceXY,self.absToleranceXY,self.absToleranceZ};
            end
        end

        function Y0 = initialConditions(self)
            % return a cell array of variables that need to be integrated
            if self.isXYOnly
                Y0 = {self.x,self.y};
            else
                Y0 = {self.x,self.y,self.z};
            end
        end

        function F = fluxAtTime(self,t,y0)
            % return a cell array of the flux of the variables being
            % integrated. You may want to call -updateIntegratorValues.
            self.updateIntegratorValues(t,y0);
            F = cell(1,self.nFluxComponents);
            if self.isXYOnly
                [F{:}] = self.wvt.variableAtPositionWithName(self.x,self.y,self.z,'u','v',interpolationMethod=self.advectionInterpolation);
            else
                [F{:}] = self.wvt.variableAtPositionWithName(self.x,self.y,self.z,'u','v','w',interpolationMethod=self.advectionInterpolation);
            end
        end

        function updateIntegratorValues(self,t,y0)
            % passes updated values of the variables being integrated.
            self.x = y0{1};
            self.y = y0{2};
            if ~self.isXYOnly
                self.z = y0{3};
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Read and write to file
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function initializeStorage(self,group)
            variables = containers.Map();

            commonKeys = {'isParticle','particleName'};
            commonVals = {uint8(1),self.name};
            attributes = containers.Map(commonKeys,commonVals);
            attributes('units') = 'unitless id number';
            attributes('particleVariableName') = 'id';
            [dim,var] = group.addDimension(strcat(self.name,'_id'),(1:self.nParticles).',attributes=attributes);
            variables('id') = var;

            % careful to create a new object each time we init
            spatialDimensionNames = self.model.wvt.spatialDimensionNames;
            for iVar=1:length(spatialDimensionNames)
                attributes = containers.Map(commonKeys,commonVals);
                attributes('units') = self.model.wvt.dimensionAnnotationWithName(spatialDimensionNames{iVar}).units;
                attributes('long_name') = strcat(self.model.wvt.dimensionAnnotationWithName(spatialDimensionNames{iVar}).description,', recorded along the particle trajectory');
                attributes('particleVariableName') = spatialDimensionNames{iVar};
                variables(spatialDimensionNames{iVar}) = group.addVariable(strcat(self.name,'_',spatialDimensionNames{iVar}),{dim.name,'t'},type="double",attributes=attributes,isComplex=false);
            end

            for iVar=1:length(self.trackedFieldNamesCell)
                varAnnotation = self.model.wvt.propertyAnnotationWithName(self.trackedFieldNamesCell{iVar});
                attributes = containers.Map(commonKeys,commonVals);
                attributes('units') = varAnnotation.units;
                attributes('long_name') = strcat(varAnnotation.description,', recorded along the particle trajectory');
                attributes('particleVariableName') = self.trackedFieldNamesCell{iVar};
                variables(self.trackedFieldNamesCell{iVar}) = group.addVariable(strcat(self.name,'_',self.trackedFieldNamesCell{iVar}),{dim.name,'t'},type="double",attributes=attributes,isComplex=false);
            end
        end

        function writeTimeStepToFile(self,group,outputIndex)
            group.variableWithName(strcat(self.name,'_x')).setValueAlongDimensionAtIndex(self.x,'t',outputIndex);
            group.variableWithName(strcat(self.name,'_y')).setValueAlongDimensionAtIndex(self.y,'t',outputIndex);
            if ~isequal(self.model.wvt.spatialDimensionNames,{'x','y'})
                group.variableWithName(strcat(self.name,'_z')).setValueAlongDimensionAtIndex(self.z,'t',outputIndex);
            end

            self.updateParticleTrackedFields();
            for iField=1:length(self.trackedFieldNamesCell)
                group.variableWithName(strcat(self.name,'_',self.trackedFieldNamesCell{iField})).setValueAlongDimensionAtIndex(self.trackedFields.(self.trackedFieldNamesCell{iField}),'t',outputIndex);
            end
        end
    end

    methods (Static)
        function os = observingSystemFromGroup(group,model,outputGroup)
            %initialize a WVObservingSystem instance from NetCDF file
            %
            % Subclasses to should override this method to enable model
            % restarts. This method works in conjunction with -writeToFile
            % to provide restart capability.
            %
            % - Topic: Initialization
            % - Declaration: os = observingSystemFromGroup(group,wvt)
            % - Parameter model: the WVModel to be used
            % - Returns os: a new instance of WVObservingSystem
            arguments
                group NetCDFGroup {mustBeNonempty}
                model WVModel {mustBeNonempty}
                outputGroup WVModelOutputGroup
            end
            % most variables will be returned with this call, but we still
            % need to fetch (x,y,z), and the tracked variables
            vars = CAAnnotatedClass.requiredPropertiesFromGroup(group);

            parentGroup = outputGroup.group;
            nPoints = parentGroup.dimensionWithName("t").nPoints;
            vars.x = parentGroup.readVariablesAtIndexAlongDimension('t',nPoints,vars.name+"_x");
            vars.y = parentGroup.readVariablesAtIndexAlongDimension('t',nPoints,vars.name+"_y");
            if ~isequal(model.wvt.spatialDimensionNames,{'x','y'})
                vars.z = parentGroup.readVariablesAtIndexAlongDimension('t',nPoints,vars.name+"_z");
            end

            options = namedargs2cell(vars);
            os = WVLagrangianParticles(model,options{:});

            % We still need to recover the last value of the tracked fields
            for i=1:length(os.trackedFieldNamesCell)
                trackedFieldName = os.name+"_"+os.trackedFieldNamesCell{i};
                os.trackedFields.(os.trackedFieldNamesCell{i}) = parentGroup.readVariablesAtIndexAlongDimension('t',nPoints,trackedFieldName);
            end
        end

        function vars = classRequiredPropertyNames()
            vars = {'name','isXYOnly','trackedFieldNames','advectionInterpolation','trackedVarInterpolation','absToleranceXY','absToleranceZ'};
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            arguments (Output)
                propertyAnnotations CAPropertyAnnotation
            end
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
            propertyAnnotations(end+1) = CAPropertyAnnotation('name','name of Lagrangian particles');
            propertyAnnotations(end+1) = CANumericProperty('isXYOnly',{}, 'bool', 'whether the advection is only applied in x-y');
            propertyAnnotations(end+1) = CAPropertyAnnotation('trackedFieldNames','tracked field names');
            propertyAnnotations(end+1) = CAPropertyAnnotation('advectionInterpolation','interpolation method for the advection scheme');
            propertyAnnotations(end+1) = CAPropertyAnnotation('trackedVarInterpolation','interpolation method for the tracked fields');
            propertyAnnotations(end+1) = CANumericProperty('absToleranceXY', {}, 'm','absolute tolerance for the adaptive integrator in x-y directions');
            propertyAnnotations(end+1) = CANumericProperty('absToleranceZ', {}, 'm','absolute tolerance for the adaptive integrator in z direction');
        end
    end
end