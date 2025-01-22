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
    properties (Access=private)
        variableAnnotationNameMap
        dimensionAnnotationNameMap
    end

    methods
        function self = PMAnnotatedClass()
            % Matlab dictionaries do not allow a subclass type; the
            % work-around is to use a cell.
            self.dimensionAnnotationNameMap = configureDictionary("string","cell");
            self.variableAnnotationNameMap = configureDictionary("string","cell");

            className = class(self);
            self.addDimensionAnnotation(feval(strcat(className,'.classDefinedDimensionAnnotations')));
            self.addVariableAnnotation(feval(strcat(className,'.classDefinedVariableAnnotations')));
        end
        
        dims = requiredDimensions(self);
        vars = requiredVariables(self);

        addDimensionAnnotation(self,dimensionAnnotation)
        val = dimensionAnnotationWithName(self,name)

        addVariableAnnotation(self,variableAnnotation)
        removeVariableAnnotation(self,variableAnnotation)
        val = variableAnnotationWithName(self,name)
        names = variableNames(self)

    end

    methods (Static,Abstract)
        dims = classRequiredDimensions()
        vars = classRequiredVariables()
    end

    methods (Static)
        dimensionAnnotations = classDefinedDimensionAnnotations()
        variableAnnotations = classDefinedVariableAnnotations()
        atc = annotatedClassFromFile(path)
    end
end