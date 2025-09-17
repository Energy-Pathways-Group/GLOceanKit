classdef WVHorizontalScalarDiffusivity < WVForcing
    % Horizontal laplacian damping with viscosity and diffusivity
    %
    % The damping is a simple horizontal Laplacian.
    %
    %
    % - Topic: Initializing
    % - Declaration: WVHorizontalScalarDiffusivity < [WVForcing](/classes/forcing/wvforcing/)
    properties
        nu
        kappa
    end

    methods
        function self = WVHorizontalScalarDiffusivity(wvt,options)
            % initialize the WVHorizontalScalarDiffusivity
            %
            % - Declaration: self = WVHorizontalScalarDiffusivity(wvt)
            % - Parameter wvt: a WVTransform instance
            % - Returns self: a WVHorizontalScalarDiffusivity instance
            arguments
                wvt WVTransform {mustBeNonempty}
                options.nu = 1e-4
                options.kappa = 1e-6
            end
            self@WVForcing(wvt,"horizontal scalar diffusivity",WVForcingType(["HydrostaticSpatial","NonhydrostaticSpatial"]));
            self.wvt = wvt;
            self.isClosure = true;
            self.nu = options.nu;
            self.kappa = options.kappa;

            % construct the damping operator
            % [K,L,~] = self.wvt.kljGrid;
            % self.F_damp = -(K.^2 +L.^2);
        end

        function [Fu, Fv, Feta] = addHydrostaticSpatialForcing(self, wvt, Fu, Fv, Feta)
            Fu = Fu + self.nu*(wvt.diffX(wvt.u,n=2) + wvt.diffY(wvt.u,n=2));
            Fv = Fv +  self.nu*(wvt.diffX(wvt.v,n=2) + wvt.diffY(wvt.v,n=2));
            Feta = Feta + self.kappa*(wvt.diffX(wvt.eta,n=2) + wvt.diffY(wvt.eta,n=2));
        end

        function [Fu, Fv, Fw, Feta] = addNonhydrostaticSpatialForcing(self, wvt, Fu, Fv, Fw, Feta)
            Fu = Fu + self.nu*(wvt.diffX(wvt.u,n=2) + wvt.diffY(wvt.u,n=2));
            Fv = Fv +  self.nu*(wvt.diffX(wvt.v,n=2) + wvt.diffY(wvt.v,n=2));
            Fw = Fw +  self.nu*(wvt.diffX(wvt.w,n=2) + wvt.diffY(wvt.w,n=2));
            Feta = Feta + self.kappa*(wvt.diffX(wvt.eta,n=2) + wvt.diffY(wvt.eta,n=2));
        end

        function force = forcingWithResolutionOfTransform(self, wvtX2)
            % Creates a forcing with the resolution of the transform
            %
            % - Declaration: forcingWithResolutionOfTransform(self, wvtX2)
            % - Parameter self: an instance of WVHorizontalScalarDiffusivity
            % - Parameter wvtX2: a WVTransform instance with doubled resolution
            % - Returns: force
            arguments
                self WVHorizontalScalarDiffusivity {mustBeNonempty}
                wvtX2 WVTransform {mustBeNonempty}
            end
            force = WVHorizontalScalarDiffusivity(wvtX2);
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
            vars = {"nu","kappa"};
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
            propertyAnnotations(end+1) = CANumericProperty('nu', {}, 'm^2 s^{-1}','viscosity');
            propertyAnnotations(end+1) = CANumericProperty('kappa', {}, 'm^2 s^{-1}','diffusivity');
        end
    end
end