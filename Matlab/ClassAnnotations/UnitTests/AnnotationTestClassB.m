classdef AnnotationTestClassB < handle & CAAnnotatedClass
    properties
        x
    end

    methods
        function self = AnnotationTestClassB(options)
            arguments
                options.x
            end
            requiredProperties = feval(strcat(class(self),'.classRequiredPropertyNames'));
            canInitializeDirectly = all(isfield(options,requiredProperties));
            if canInitializeDirectly
                for iVar = 1:length(requiredProperties)
                    name = requiredProperties{iVar};
                    self.(name) = options.(name);
                end
            else
                self.x = linspace(-100,0,101);
            end
        end
    end

    methods (Static)

        function vars = classRequiredPropertyNames()
            vars = {'x'};
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            arguments (Output)
                propertyAnnotations CAPropertyAnnotation
            end
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
            propertyAnnotations(end+1) = CADimensionProperty('x', 'm', 'x coordinate');
        end
        
        function atc = annotatedTestClassFromFile(path)
            atc = CAAnnotatedClass.annotatedClassFromFile(path);
        end

    end
end