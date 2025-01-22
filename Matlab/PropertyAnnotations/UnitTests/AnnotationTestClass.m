classdef AnnotationTestClass < handle & PMAnnotatedClass
    properties
        z
        f
        myVar
    end

    methods
        function self = AnnotationTestClass(options)
            arguments
                options.z
                options.f
                options.myVar
            end
            className = class(self);
            requiredVariables = union(feval(strcat(className,'.classRequiredDimensions')),feval(strcat(className,'.classRequiredVariables')));
            canInitializeDirectly = all(isfield(options,requiredVariables));
            if canInitializeDirectly
                for iVar = 1:length(requiredVariables)
                    name = requiredVariables{iVar};
                    self.(name) = options.(name);
                end
            else
                self.z = linspace(-100,0,101);
                self.f = @(z) exp(z/20);
                self.myVar = self.z.^2;
            end
        end
    end

    methods (Static)
        function dims = classRequiredDimensions()
            dims = {'z'};
        end

        function vars = classRequiredVariables()
            vars = {'f','myVar'};
        end

        function dimensions = classDefinedDimensionAnnotations()
            dimensions = PMDimensionAnnotation.empty(0,0);

            dimensions(end+1) = PMDimensionAnnotation('z', 'm', 'z coordinate');
            dimensions(end).attributes('standard_name') = 'height_above_mean_sea_level';
            dimensions(end).attributes('positive') = 'up';
            dimensions(end).attributes('axis') = 'Z';
        end

        function variableAnnotations = classDefinedVariableAnnotations()
            variableAnnotations = PMVariableAnnotation.empty(0,0);
            variableAnnotations(end+1) = PMVariableAnnotation('myVar',{'z'},'', 'A variable quadratic in z.');
            variableAnnotations(end+1) = PMVariableAnnotation('f',{},'', 'A function handle that does something!');
        end
        
        function atc = annotatedTestClassFromFile(path)
            atc = PMAnnotatedClass.annotatedClassFromFile(path);
        end

    end
end