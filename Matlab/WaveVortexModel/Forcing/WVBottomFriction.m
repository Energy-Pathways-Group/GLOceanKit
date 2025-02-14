classdef WVBottomFriction < WVForcing
    % Bottom friction
    %
    % Applies linear bottom friction to the flow, i.e., $$\frac{du}{dt} = -r*u$$.
    %
    % - Topic: Initializing
    % - Declaration: WVBottomFriction < [WVForcing](/classes/wvforcing/)
    properties
        r
        RVA0
    end

    methods
        function self = WVBottomFriction(wvt,options)
            % initialize the WVNonlinearFlux nonlinear flux
            %
            % - Declaration: self = WVBottomFriction(wvt,options)
            % - Parameter wvt: a WVTransform instance
            % - Parameter r: (optional) linear bottom friction, try 1/(200*86400)
            % - Returns frictionalForce: a WVBottomFriction instance
            arguments
                wvt WVTransform {mustBeNonempty}
                options.r (1,1) double {mustBeNonnegative} = 0 % linear bottom friction, try 1/(200*86400) https://www.nemo-ocean.eu/doc/node70.html
            end
            self@WVForcing(wvt,"bottom friction",WVForcingType(["HydrostaticSpatial" "NonhydrostaticSpatial" "PVSpatial"]));
            self.r = options.r;
            self.RVA0 = - wvt.g * (wvt.Kh .* wvt.Kh) / wvt.f;
        end

        function [Fu, Fv, Feta] = addHydrostaticSpatialForcing(self, wvt, Fu, Fv, Feta)
            Fu(:,:,1) = Fu(:,:,1) - self.r*wvt.u(:,:,1);
            Fv(:,:,1) = Fv(:,:,1) - self.r*wvt.v(:,:,1);
        end

        function [Fu, Fv, Fw, Feta] = addNonhydrostaticSpatialForcing(self, wvt, Fu, Fv, Fw, Feta)
            Fu(:,:,1) = Fu(:,:,1) - self.r*wvt.u(:,:,1);
            Fv(:,:,1) = Fv(:,:,1) - self.r*wvt.v(:,:,1);
        end

        function Fpv = addPotentialVorticitySpatialForcing(self, wvt, Fpv)
            mask = zeros(wvt.spatialMatrixSize);
            mask(:,:,1) = 1;
            Fpv = Fpv - self.r * mask.*wvt.transformToSpatialDomainWithF(A0=self.RVA0 .* wvt.A0);
        end

        function F0 = addPotentialVorticitySpectralForcing(self, wvt, F0)
            F0 = F0 + self.damp .* wvt.A0;
        end

        function force = forcingWithResolutionOfTransform(self,wvtX2)
            force = WVBottomFriction(wvtX2,r=self.r);
        end
    end
    methods (Static)
        function vars = classRequiredPropertyNames()
            vars = {'r'};
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            arguments (Output)
                propertyAnnotations CAPropertyAnnotation
            end
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
            propertyAnnotations(end+1) = CANumericProperty('r', {}, 's^{-1}','bottom friction');
        end
    end

end