classdef WVAnnotatedClass < handle
    %WVAnnotatedClass A class with metadata
    %
    % An `annotated class' has metadata describing its methods and
    % properties. This metadata is used to read/write from file, and to
    % create website documentation.
    %
    % Annotations should specify the netcdf group they belong to.
    % Annotations should specify whether the type is a function_handle
    %
    % A subclass must do two things:
    % 1) implement the required methods (perhaps by simply calling the
    % static methods) and,
    % 2) initialize the superclass during init, self@WVAnnotatedClass()
    % 
    properties %(Access=private)
        variableAnnotationNameMap
        timeDependentVariablesNameMap
        propertyAnnotationNameMap
        dimensionAnnotationNameMap
    end

    methods (Abstract)
        dims = requiredDimensions(self);
        vars = requiredVariables(self);
        dimAnnotations = dimensionAnnotations(self);
        propAnnotations = propertyAnnotations(self);
        varAnnotations = variableAnnotations(self);
    end

    methods
        function self = WVAnnotatedClass()
            self.dimensionAnnotationNameMap = configureDictionary("string","WVDimensionAnnotation");
            self.propertyAnnotationNameMap = configureDictionary("string","WVPropertyAnnotation"); 
            self.variableAnnotationNameMap = configureDictionary("string","WVVariableAnnotation");
            self.timeDependentVariablesNameMap = configureDictionary("string","WVVariableAnnotation");

            self.addDimensionAnnotations(self.dimensionAnnotations);
            self.addPropertyAnnotations(self.propertyAnnotations);
            self.addVariableAnnotations(self.variableAnnotations);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Metadata
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        addDimensionAnnotations(self,dimensionAnnotation)
        val = dimensionAnnotationWithName(self,name)

        addPropertyAnnotations(self,propertyAnnotation)
        val = propertyAnnotationWithName(self,name)

        addVariableAnnotations(self,variableAnnotation)
        removeVariableAnnotations(self,variableAnnotation)
        val = variableAnnotationWithName(self,name)
        names = variableNames(self)

    end

    methods (Static)
        function matFilePath = matlabSidecarPathForNetCDFPath(path)
            [filepath,name,~] = fileparts(path);
            if isempty(filepath)
                matFilePath = sprintf('%s.mat',name);
            else
                matFilePath = sprintf('%s/%s.mat',filepath,name);
            end
        end

        function dims = classRequiredDimensions()
            dims = {};
        end

        function vars = classRequiredVariables()
            vars = {};
        end

        function dimensions = classDimensionAnnotations()
            % return array of WVDimensionAnnotation to annotate the
            % dimensions
            %
            % This function returns annotations for all dimensions of the
            % WVStratification class.
            %
            % - Topic: Developer
            % - Declaration: dimensionAnnotations = WVStratification.dimensionAnnotationsForStratifiedFlow()
            % - Returns dimensionAnnotations: array of WVDimensionAnnotation instances
            arguments (Output)
                dimensions WVDimensionAnnotation
            end
            dimensions = WVDimensionAnnotation.empty(0,0);
        end

        function propertyAnnotations = classPropertyAnnotations()
            % return array of WVPropertyAnnotation initialized by default
            %
            % This function returns annotations for all properties of the
            % WVStratification class.
            %
            % - Topic: Developer
            % - Declaration: propertyAnnotations = WVStratification.propertyAnnotationsForStratifiedFlow()
            % - Returns propertyAnnotations: array of WVPropertyAnnotation instances
            propertyAnnotations = WVPropertyAnnotation.empty(0,0);
        end

        function methodAnnotations = classMethodAnnotations()
            % return array of WVAnnotations to annotate the methods
            %
            % This function returns annotations for all methods of the
            % WVStratification class.
            %
            % - Topic: Developer
            % - Declaration: methodAnnotations = WVStratification.methodAnnotationsForStratifiedFlow()
            % - Returns methodAnnotations: array of WVAnnotations instances
            methodAnnotations = WVAnnotation.empty(0,0);
        end
    end
end