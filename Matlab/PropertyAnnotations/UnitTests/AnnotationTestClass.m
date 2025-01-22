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
            requiredProperties = union(feval(strcat(className,'.classRequiredDimensions')),feval(strcat(className,'.classRequiredProperties')));
            canInitializeDirectly = all(isfield(options,requiredProperties));
            if canInitializeDirectly
                for iVar = 1:length(requiredProperties)
                    name = requiredProperties{iVar};
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

        function vars = classRequiredProperties()
            vars = {'f','myVar'};
        end

        function dimensions = classDefinedDimensionAnnotations()
            dimensions = PMDimensionAnnotation.empty(0,0);

            dimensions(end+1) = PMDimensionAnnotation('z', 'm', 'z coordinate');
            dimensions(end).attributes('standard_name') = 'height_above_mean_sea_level';
            dimensions(end).attributes('positive') = 'up';
            dimensions(end).attributes('axis') = 'Z';
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = PMPropertyAnnotation.empty(0,0);
            propertyAnnotations(end+1) = PMPropertyAnnotation('myVar',{'z'},'', 'A variable quadratic in z.');
            propertyAnnotations(end+1) = PMPropertyAnnotation('f',{},'', 'A function handle that does something!');
        end
        
        function atc = annotatedTestClassFromFile(path)
            atc = PMAnnotatedClass.annotatedClassFromFile(path);
        end

    end
end