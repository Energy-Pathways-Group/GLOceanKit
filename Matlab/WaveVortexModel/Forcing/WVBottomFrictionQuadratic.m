classdef WVBottomFrictionQuadratic < WVForcing
    % Quadratic bottom friction
    %
    % Applies quadratic bottom friction to the flow, i.e.,
    % $$\frac{du}{dt} = -(Cd/dz)*|u|*u$$. Cd is unitless, and dz is
    % (approximately) the size of the grid spacing at the bottom boundary.
    %
    % To compare with linear bottom friction where $$\frac{du}{dt} =
    % -r*u$$, note that $$-(Lz/dz)*r = -(Cd/dz)*|u|$$ and you will find a
    % characteristic velocity $$|u|$$ of about 11.5 cm/s for Cd=0.002 and
    % r=1/(200 days). If Cd=0.001, then the damping time scale has to
    % double to 1/(400 days) for and equivalent characteristic velocity.
    %
    % For barotropic QG we want the scaling to work out similarly, but now
    % have $$-r = -(Cd/D)*|u|$$ where $$D$$ works out to be Lz. Thus, we
    % will scale the barotropic QG quadratic bottom drag by 4000 m, to match
    % typical oceanic scales.
    % 
    %
    % - Topic: Initializing
    % - Declaration: WVBottomFrictionQuadratic < [WVForcing](/classes/wvforcing/)
    properties
        Cd % non-dimensional
        cd % units of inverse length
    end

    methods
        function self = WVBottomFrictionQuadratic(wvt,options)
            % initialize the WVNonlinearFlux nonlinear flux
            %
            % - Declaration: self = WVBottomFrictionQuadratic(wvt,options)
            % - Parameter wvt: a WVTransform instance
            % - Parameter Cd: (optional) non-dimensional quadratic damping coefficient. Default is 0.001
            % - Returns frictionalForce: a WVBottomFriction instance
            arguments
                wvt WVTransform {mustBeNonempty}
                options.Cd (1,1) double {mustBeNonnegative} = 1e-3 % https://www.nemo-ocean.eu/doc/node70.html
            end
            self@WVForcing(wvt,"quadratic bottom friction",WVForcingType(["HydrostaticSpatial" "NonhydrostaticSpatial" "PVSpatial"]));
            self.Cd = options.Cd;
            
            if ~isa(self.wvt,"WVGeometryDoublyPeriodicBarotropic")
                self.cd = self.Cd/wvt.z_int(1);
            else
                % scaled by Lz, a typical ocean depth
                self.cd = self.Cd/4000;
            end
        end

        function [Fu, Fv, Feta] = addHydrostaticSpatialForcing(self, wvt, Fu, Fv, Feta)
            ub = wvt.u(:,:,1);
            vb = wvt.v(:,:,1);
            cb = sqrt(ub.^2 + vb.^2);
            Fu(:,:,1) = Fu(:,:,1) - self.cd*ub.*cb;
            Fv(:,:,1) = Fv(:,:,1) - self.cd*vb.*cb;
        end

        function [Fu, Fv, Fw, Feta] = addNonhydrostaticSpatialForcing(self, wvt, Fu, Fv, Fw, Feta)
            ub = wvt.u(:,:,1);
            vb = wvt.v(:,:,1);
            cb = sqrt(ub.^2 + vb.^2);
            Fu(:,:,1) = Fu(:,:,1) - self.cd*ub.*cb;
            Fv(:,:,1) = Fv(:,:,1) - self.cd*vb.*cb;
        end

        function Fpv = addPotentialVorticitySpatialForcing(self, wvt, Fpv)
            u_b = wvt.u(:,:,1);
            v_b = wvt.v(:,:,1);
            uv_mag = sqrt(u_b.^2 + v_b.^2);
            Fpv(:,:,1) = Fpv(:,:,1) - self.cd * (wvt.diffX(uv_mag.*v_b) - wvt.diffY(uv_mag.*u_b));
        end

        function force = forcingWithResolutionOfTransform(self,wvtX2)
            force = WVBottomFrictionQuadratic(wvtX2,Cd=self.Cd);
        end
    end
    methods (Static)
        function vars = classRequiredPropertyNames()
            vars = {'Cd'};
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            arguments (Output)
                propertyAnnotations CAPropertyAnnotation
            end
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
            propertyAnnotations(end+1) = CANumericProperty('Cd', {}, '','non-dimensional quadratic drag coefficient');
        end
    end

end