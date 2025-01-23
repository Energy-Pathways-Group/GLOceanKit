classdef PMAnnotatedClass < handle
    %PMAnnotatedClass A class with metadata describing its properties
    %
    % An `annotated class' has metadata describing its properties as either
    % dimensions or properties (which may have dimensions). This metadata
    % is used to read/write NetCDF files, and to create website
    % documentation.
    %
    % Annotations should specify the netcdf group they belong to.
    % Annotations should specify whether the type is a function_handle
    %
    % A subclass must do two things:
    % 1) implement the required methods (perhaps by simply calling the
    % static methods) and,
    % 2) initialize the superclass during init, self@PMAnnotatedClass()
    % 
    properties (Access=private)
        propertyAnnotationNameMap
        dimensionAnnotationNameMap
    end

    methods
        function self = PMAnnotatedClass()
            % Matlab dictionaries do not allow a subclass type; the
            % work-around is to use a cell.
            self.dimensionAnnotationNameMap = configureDictionary("string","cell");
            self.propertyAnnotationNameMap = configureDictionary("string","cell");

            className = class(self);
            self.addDimensionAnnotation(feval(strcat(className,'.classDefinedDimensionAnnotations')));
            self.addPropertyAnnotation(feval(strcat(className,'.classDefinedPropertyAnnotations')));
        end
        
        dims = requiredDimensions(self);
        vars = requiredProperties(self);

        addDimensionAnnotation(self,dimensionAnnotation)
        val = dimensionAnnotationWithName(self,name)
        function dimAnnotations = dimensionAnnotations(self)
            arguments (Input)
                self PMAnnotatedClass
            end
            arguments (Output)
                dimAnnotations PMDimensionAnnotation
            end
            dimAnnotations = [self.dimensionAnnotationNameMap{self.dimensionAnnotationNameMap.keys}];
        end

        addPropertyAnnotation(self,variableAnnotation)
        removePropertyAnnotation(self,variableAnnotation)
        val = propertyAnnotationWithName(self,name)
        names = propertyNames(self)
        function propAnnotations = propertyAnnotations(self)
            arguments (Input)
                self PMAnnotatedClass
            end
            arguments (Output)
                propAnnotations PMPropertyAnnotation
            end
            propAnnotations = [self.propertyAnnotationNameMap{self.propertyAnnotationNameMap.keys}];
        end

        ncfile = writeToFile(self,path,variables,options)
        
    end

    methods (Static,Abstract)
        dims = classRequiredDimensions()
        vars = classRequiredProperties()
    end

    methods (Static)
        dimensionAnnotations = classDefinedDimensionAnnotations()
        propertyAnnotations = classDefinedPropertyAnnotations()
        atc = annotatedClassFromFile(path)

        function var = requiredPropertiesFromGroup(group,options)
            arguments (Input)
                group NetCDFGroup
                options.className string
            end
            arguments (Output)
                var struct
            end
            if isfield(options,'className')
                className = options.className;
            else
                if isKey(group.attributes,'AnnotatedClass')
                    className = group.attributes('AnnotatedClass');
                else
                    error('Unable to find the attribute AnnotatedClass');
                end
            end
            requiredProperties = union(feval(strcat(className,'.classRequiredDimensions')),feval(strcat(className,'.classRequiredProperties')));
            var = PMAnnotatedClass.variablesFromGroup(group,requiredProperties);
        end

        function var = variablesFromGroup(group,variables)
            arguments (Input)
                group NetCDFGroup
                variables {mustBeA(variables,"cell")} = {}
            end
            arguments (Output)
                var struct
            end
            for iVar = 1:length(variables)
                name = variables{iVar};
                var.(name) = group.readVariables(name);
            end
        end

        function [canInit, errorString] = canInitializeDirectlyFromGroup(group,requiredDimensions,requiredProperties)
            canInit = false;
            errorString = "";
            hasDim = group.hasVariableWithName(requiredDimensions{:});
            if ~all(hasDim)
                errorString = "The NetCDF file is missing required dimensions: " + join(string(requiredDimensions(~hasDim)),', ');
                return;
            end

            hasVar = group.hasVariableWithName(requiredProperties{:});
            if ~all(hasVar)
                errorString = "The NetCDF file is missing required variables: " + join(string(requiredProperties(~hasVar)),', ');
                return;
            end

            canInit = true;
        end

        function ncfile = writeToPath(ac,path,options)
            arguments (Input)
                ac PMAnnotatedClass {mustBeNonempty}
                path char {mustBeNonempty}
                options.shouldOverwriteExisting double {mustBeMember(options.shouldOverwriteExisting,[0 1])} = 0
                options.dims {mustBeA(options.dims,"cell")} = {}
                options.properties {mustBeA(options.properties,"cell")} = {}
                options.dimAnnotations PMDimensionAnnotation = PMDimensionAnnotation.empty(0,0);
                options.propAnnotations PMPropertyAnnotation = PMPropertyAnnotation.empty(0,0);
                options.attributes = configureDictionary("string","string")
            end
            arguments (Output)
                ncfile NetCDFFile
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
            ncfile = NetCDFFile(path);
            PMAnnotatedClass.writeToGroup(ac,ncfile,dims=options.dims,properties=options.properties,dimAnnotations=options.dimAnnotations,propAnnotations=options.propAnnotations,attributes=options.attributes);
        end

        function writeToGroup(ac,group,options)
            arguments
                ac PMAnnotatedClass
                group NetCDFGroup
                options.dims {mustBeA(options.dims,"cell")} = {}
                options.properties {mustBeA(options.properties,"cell")} = {}
                options.dimAnnotations PMDimensionAnnotation = PMDimensionAnnotation.empty(0,0);
                options.propAnnotations PMPropertyAnnotation = PMPropertyAnnotation.empty(0,0);
                options.attributes = configureDictionary("string","string")
            end

            group.addAttribute('AnnotatedClass',class(ac));

            attributeNames = keys(options.attributes);
            for iKey=1:length(attributeNames)
                group.addAttribute(attributeNames{iKey},options.attributes(attributeNames{iKey}));
            end

            dimensionAnnotationDictionary = dictionary(string({options.dimAnnotations.name}),num2cell(options.dimAnnotations));
            for iDim=1:length(options.dims)
                dimAnnotation = dimensionAnnotationDictionary{options.dims{iDim}};
                dimAnnotation.attributes('units') = dimAnnotation.units;
                dimAnnotation.attributes('long_name') = dimAnnotation.description;
                group.addDimension(dimAnnotation.name,ac.(dimAnnotation.name),attributes=dimAnnotation.attributes);
            end

            propertyAnnotationDictionary = dictionary(string({options.propAnnotations.name}),num2cell(options.propAnnotations));
            for iVar=1:length(options.properties)
                % check for function_handle, and add group
                propAnnotation = propertyAnnotationDictionary{options.properties{iVar}};
                if ~isempty(propAnnotation.units)
                    propAnnotation.attributes('units') = propAnnotation.units;
                end
                if ~isempty(propAnnotation.description)
                    propAnnotation.attributes('long_name') = propAnnotation.description;
                end
                if isa(ac.(propAnnotation.name),'function_handle')
                    group.addFunctionHandle(propAnnotation.name,ac.(propAnnotation.name),attributes=propAnnotation.attributes);
                else
                    group.addVariable(propAnnotation.name,propAnnotation.dimensions,ac.(propAnnotation.name),isComplex=propAnnotation.isComplex,attributes=propAnnotation.attributes);
                end
            end
        end
    end
end