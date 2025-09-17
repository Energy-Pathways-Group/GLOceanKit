classdef WVVerticalDiffusivity < WVForcing
    % Vertical diffusivty
    %
    %
    % - Topic: Initializing
    % - Declaration: WVVerticalDiffusivity < [WVForcing](/classes/forcing/wvforcing/)
    properties
        kappa_z
        shouldForceMeanDensityAnomaly
        dLnN2 = 0
    end

    methods
        function self = WVVerticalDiffusivity(wvt,options)
            % initialize the WVNonlinearFlux nonlinear flux
            %
            % - Declaration: nlFlux = WVAdaptiveViscosity(wvt)
            % - Parameter wvt: a WVTransform instance
            % - Returns self: a WVAdaptiveViscosity instance
            arguments
                wvt WVTransform {mustBeNonempty}
                options.kappa_z double = 1e-5
                options.shouldForceMeanDensityAnomaly = true;
            end
            self@WVForcing(wvt,"vertical diffusivity",WVForcingType(["HydrostaticSpatial","NonhydrostaticSpatial","PVSpatial"]));
            self.wvt = wvt;
            self.kappa_z = options.kappa_z;
            self.shouldForceMeanDensityAnomaly = options.shouldForceMeanDensityAnomaly;
            if isa(wvt,'WVStratificationVariable') && self.shouldForceMeanDensityAnomaly
                self.dLnN2 = shiftdim(wvt.dLnN2,-2);
            end
        end

        function [Fu, Fv, Feta] = addHydrostaticSpatialForcing(self, wvt, Fu, Fv, Feta)
            Feta = Feta + self.kappa_z * (wvt.diffZG(wvt.eta,n=2) - self.dLnN2);
        end

        function [Fu, Fv, Fw, Feta] = addNonhydrostaticSpatialForcing(self, wvt, Fu, Fv, Fw, Feta)
            Feta = Feta + self.kappa_z * (wvt.diffZG(wvt.eta,n=2) - self.dLnN2);
        end

        function Fpv = addPotentialVorticitySpatialForcing(self, wvt, Fpv)
            % Fpv = Fpv - wvt.f * self.kappa_z * (wvt.diffZG(wvt.eta,3) - wvt.diffZG(self.dLnN2));
            Fpv = Fpv - wvt.f * self.kappa_z * (wvt.diffZG(wvt.eta,n=3));
            % I believe this is incorrect because it excludes
            % the MDA
        end

        function force = forcingWithResolutionOfTransform(self, wvtX2)
            % Creates a forcing with the resolution of the transform
            %
            % - Declaration: forcingWithResolutionOfTransform(self, wvtX2)
            % - Parameter self: an instance of WVAdaptiveViscosity
            % - Parameter wvtX2: a WVTransform instance with doubled resolution
            % - Returns: force
            arguments
                self WVVerticalDiffusivity {mustBeNonempty}
                wvtX2 WVTransform {mustBeNonempty}
            end
            force = WVVerticalDiffusivity(wvtX2);
        end
    end
    methods (Static)
        function vars = classRequiredPropertyNames()
            % Returns the required property names for the class
            %
            % - Declaration: classRequiredPropertyNames()
            % - Returns: vars
            arguments
            end
            vars = {"kappa_z","shouldForceMeanDensityAnomaly"};
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            % Returns the defined property annotations for the class
            %
            % - Declaration: classDefinedPropertyAnnotations()
            % - Returns: propertyAnnotations
            arguments (Output)
                propertyAnnotations CAPropertyAnnotation
            end
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
            propertyAnnotations(end+1) = CANumericProperty('kappa_z', {}, 'm^2 s^{-1}','vertical diffusivity');
            propertyAnnotations(end+1) = CANumericProperty('shouldForceMeanDensityAnomaly',{},'bool', 'whether the vertical diffusivity is applied to the mean density anomaly');
        end
    end
end