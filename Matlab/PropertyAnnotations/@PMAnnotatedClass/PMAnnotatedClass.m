classdef PMAnnotatedClass < handle
    %PMAnnotatedClass A class with metadata describing its properties
    %
    % An `annotated class' has metadata describing its properties as either
    % dimensions or variables. This metadata is used to read/write NetCDF
    % files, and to create website documentation.
    %
    % Annotations should specify the netcdf group they belong to.
    % Annotations should specify whether the type is a function_handle
    %
    % A subclass must do two things:
    % 1) implement the required methods (perhaps by simply calling the
    % static methods) and,
    % 2) initialize the superclass during init, self@PMAnnotatedClass()
    % 
    properties %(Access=private)
        variableAnnotationNameMap
        dimensionAnnotationNameMap
    end

    methods (Abstract)
        dims = requiredDimensions(self);
        vars = requiredVariables(self);
        dimAnnotations = dimensionAnnotations(self);
        varAnnotations = variableAnnotations(self);
    end

    methods
        function self = PMAnnotatedClass()
            % Matlab dictionaries do not allow a subclass type; the
            % work-around is to use a cell.
            self.dimensionAnnotationNameMap = configureDictionary("string","cell");
            self.variableAnnotationNameMap = configureDictionary("string","cell");

            self.addDimensionAnnotation(self.dimensionAnnotations);
            self.addVariableAnnotation(self.variableAnnotations);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Metadata
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        addDimensionAnnotation(self,dimensionAnnotation)
        val = dimensionAnnotationWithName(self,name)

        addVariableAnnotation(self,variableAnnotation)
        removeVariableAnnotation(self,variableAnnotation)
        val = variableAnnotationWithName(self,name)
        names = variableNames(self)

    end

    methods (Static)
        function dims = classRequiredDimensions()
            dims = {};
        end

        function vars = classRequiredVariables()
            vars = {};
        end

        function dimensions = classDimensionAnnotations()
            % return array of PMDimensionAnnotation to annotate the
            % dimensions
            %
            % This function returns annotations for all dimensions of the
            % PMAnnotatedClass class.
            %
            % - Topic: Developer
            % - Declaration: dimensionAnnotations = PMAnnotatedClass.classDimensionAnnotations()
            % - Returns dimensionAnnotations: array of PMDimensionAnnotation instances
            arguments (Output)
                dimensions PMDimensionAnnotation
            end
            dimensions = PMDimensionAnnotation.empty(0,0);
        end

        function variableAnnotations = classVariableAnnotations()
            % return array of WVPropertyAnnotation initialized by default
            %
            % This function returns annotations for all properties of the
            % PMAnnotatedClass class.
            %
            % - Topic: Developer
            % - Declaration: propertyAnnotations = PMAnnotatedClass.classVariableAnnotations()
            % - Returns propertyAnnotations: array of WVPropertyAnnotation instances
            variableAnnotations = PMVariableAnnotation.empty(0,0);
        end

        function atc = annotatedClassFromFile(path)
            ncfile = NetCDFFile(path);
            if isKey(ncfile.attributes,'PMAnnotatedClass')
                className = ncfile.attributes('PMAnnotatedClass');
                requiredVariables = union(feval(strcat(className,'.classRequiredDimensions')),feval(strcat(className,'.classRequiredVariables')));
                for iVar = 1:length(requiredVariables)
                    name = requiredVariables{iVar};
                    var.(name) = ncfile.readVariables(name);
                end
                varCell = namedargs2cell(var);
                atc = feval(className,varCell{:});
            else
                error('Unable to find the attribute PMAnnotatedClass');
            end
        end
    end
end