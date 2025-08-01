classdef AnnotationTestClass < handle & CAAnnotatedClass
    properties
        z
        f
        myVar
        myObjs
    end

    methods
        function self = AnnotationTestClass(options)
            arguments
                options.z
                options.f
                options.myVar
                options.myObjs
            end
            requiredProperties = feval(strcat(class(self),'.classRequiredPropertyNames'));
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
                self.myObjs = AnnotationTestClassB();
            end
        end
    end

    methods (Static)

        function vars = classRequiredPropertyNames()
            vars = {'z','f','myVar','myObjs'};
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            arguments (Output)
                propertyAnnotations CAPropertyAnnotation
            end
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);

            propertyAnnotations(end+1) = CADimensionProperty('z', 'm', 'z coordinate');
            propertyAnnotations(end).attributes('standard_name') = 'height_above_mean_sea_level';
            propertyAnnotations(end).attributes('positive') = 'up';
            propertyAnnotations(end).attributes('axis') = 'Z';

            propertyAnnotations(end+1) = CANumericProperty('myVar',{'z'},'', 'A variable quadratic in z.');
            propertyAnnotations(end+1) = CAFunctionProperty('f', 'A function handle that does something!');
            propertyAnnotations(end+1) = CAObjectProperty('myObjs', 'A bunch of objects');
        end
        
        function atc = annotatedTestClassFromFile(path)
            atc = CAAnnotatedClass.annotatedClassFromFile(path);
        end

    end
end