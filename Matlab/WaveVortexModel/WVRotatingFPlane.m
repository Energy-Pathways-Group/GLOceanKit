classdef WVRotatingFPlane < handle
    properties
        planetaryRadius
        rotationRate
        latitude
        g
    end
    properties (Dependent, SetAccess=private)
        inertialPeriod, f, beta
    end

    methods
        function self = WVRotatingFPlane(rotatingOptions)
            arguments
                rotatingOptions.rotationRate (1,1) double = 7.2921E-5
                rotatingOptions.planetaryRadius (1,1) double = 6.371e6
                rotatingOptions.latitude (1,1) double = 33
                rotatingOptions.g (1,1) double = 9.81
            end
            self.planetaryRadius = rotatingOptions.planetaryRadius;
            self.rotationRate = rotatingOptions.rotationRate;
            self.latitude = rotatingOptions.latitude;
            self.g = rotatingOptions.g;
        end
        function value = get.inertialPeriod(self)
            value = (2*pi/(2 * self.rotationRate * sin( self.latitude*pi/180 )));
        end

        function value = get.f(self)
            value = 2 * self.rotationRate * sin( self.latitude*pi/180 );
        end

        function value = get.beta(self)
            value = 2 * self.rotationRate * cos( self.latitude*pi/180 ) / self.planetaryRadius;
        end
    end

    methods (Static)
        function requiredPropertyNames = namesOfRequiredPropertiesForRotatingFPlane()
            requiredPropertyNames = {'planetaryRadius','rotationRate','g','latitude'};
        end

        function propertyAnnotations = propertyAnnotationsForRotatingFPlane()
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
            propertyAnnotations(end+1) = CANumericProperty('latitude',{},'degrees_north', 'central latitude of the simulation', detailedDescription='- topic: Domain Attributes');
            propertyAnnotations(end).attributes('standard_name') = 'latitude';
            propertyAnnotations(end+1) = CANumericProperty('rotationRate',{},'rad/s', 'rotation rate of the planetary body', detailedDescription='- topic: Domain Attributes');
            propertyAnnotations(end+1) = CANumericProperty('planetaryRadius',{},'m', 'radius of the planetary body', detailedDescription='- topic: Domain Attributes');
            propertyAnnotations(end+1) = CANumericProperty('f',{},'rad/s', 'Coriolis parameter', detailedDescription='- topic: Domain Attributes');
            propertyAnnotations(end+1) = CANumericProperty('inertialPeriod',{},'s', 'inertial period');
            propertyAnnotations(end+1) = CANumericProperty('g',{},'m s^{-2}', 'gravity of Earth', detailedDescription='- topic: Domain Attributes');
        end

        function rotatingOptions = requiredPropertiesForRotatingFPlaneFromGroup(group)
            arguments (Input)
                group NetCDFGroup {mustBeNonempty}
            end
            arguments (Output)
                rotatingOptions
            end
            requiredProperties = WVRotatingFPlane.namesOfRequiredPropertiesForRotatingFPlane;
            [canInit, errorString] = CAAnnotatedClass.canInitializeDirectlyFromGroup(group,requiredProperties);
            if ~canInit
                error(errorString);
            end

            vars = CAAnnotatedClass.propertyValuesFromGroup(group,requiredProperties);
            rotatingOptions = namedargs2cell(vars);
        end
    end
end