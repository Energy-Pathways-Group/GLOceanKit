classdef WVBottomFrictionLinear < WVForcing
    % Linear bottom friction
    %
    % Applies linear bottom friction to the flow, i.e., $$\frac{du}{dt} = -r*u$$.
    %
    % The linear bottom friction is scaled such that we actually apply,
    % $$\frac{du}{dt} = -(Lz/dz)*r*u$$ and the volume integrated effect of
    % friction remains the same regardless of resolution.
    %
    % To compare with quadratic bottom friction where $$\frac{du}{dt} =
    % -(Cd/dz)*|u|*u$$, note that $$-(Lz/dz)*r = -(Cd/dz)*|u|$$ and you
    % will find a characteristic velocity $$|u|$$ of about 10 cm/s for
    % Cd=0.002.
    %
    % - Topic: Initializing
    % - Declaration: WVBottomFriction < [WVForcing](/classes/wvforcing/)
    properties
        r
        r_scaled
        % RVA0
    end

    methods
        function self = WVBottomFrictionLinear(wvt,options)
            % initialize the WVNonlinearFlux nonlinear flux
            %
            % - Declaration: self = WVBottomFriction(wvt,options)
            % - Parameter wvt: a WVTransform instance
            % - Parameter r: (optional) linear bottom friction, try 1/(200*86400)
            % - Returns frictionalForce: a WVBottomFriction instance
            arguments
                wvt WVTransform {mustBeNonempty}
                options.r (1,1) double {mustBeNonnegative} = 1/(200*86400) % linear bottom friction, try 1/(200*86400) https://www.nemo-ocean.eu/doc/node70.html
            end
            self@WVForcing(wvt,"linear bottom friction",WVForcingType(["HydrostaticSpatial" "NonhydrostaticSpatial" "PVSpatial"]));
            self.r = options.r;
            if ~isa(self.wvt,"WVGeometryDoublyPeriodicBarotropic")
                self.r_scaled = self.r * wvt.Lz / wvt.z_int(1);
            else
                self.r_scaled = self.r;
            end
        end

        function [Fu, Fv, Feta] = addHydrostaticSpatialForcing(self, wvt, Fu, Fv, Feta)
            Fu(:,:,1) = Fu(:,:,1) - self.r_scaled*wvt.u(:,:,1);
            Fv(:,:,1) = Fv(:,:,1) - self.r_scaled*wvt.v(:,:,1);
        end

        function [Fu, Fv, Fw, Feta] = addNonhydrostaticSpatialForcing(self, wvt, Fu, Fv, Fw, Feta)
            Fu(:,:,1) = Fu(:,:,1) - self.r_scaled*wvt.u(:,:,1);
            Fv(:,:,1) = Fv(:,:,1) - self.r_scaled*wvt.v(:,:,1);
        end

        function Fpv = addPotentialVorticitySpatialForcing(self, wvt, Fpv)
            Fpv(:,:,1) = Fpv(:,:,1) - self.r_scaled * wvt.zeta_z(:,:,1);
        end

        function force = forcingWithResolutionOfTransform(self,wvtX2)
            force = WVBottomFrictionLinear(wvtX2,r=self.r);
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