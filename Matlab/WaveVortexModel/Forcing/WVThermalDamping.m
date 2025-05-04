classdef WVThermalDamping < WVForcing
    % Thermal damping
    %
    % Applies thermal damping to the flow, i.e., $$\frac{dq}{dt} = \alpha \lambda^2 \psi$$.
    %
    % This is as defined in Scott and Dritschel, but it can be shown that
    % it is basically just a vertical diffusivity.
    %
    % - Topic: Initializing
    % - Declaration: WVBottomFriction < [WVForcing](/classes/wvforcing/)
    properties
        alpha
        alpha_scaled
        % RVA0
    end

    methods
        function self = WVThermalDamping(wvt,options)
            % initialize the WVNonlinearFlux nonlinear flux
            %
            % - Declaration: self = WVBottomFriction(wvt,options)
            % - Parameter wvt: a WVTransform instance
            % - Parameter r: (optional) linear bottom friction, try 1/(200*86400)
            % - Returns frictionalForce: a WVBottomFriction instance
            arguments
                wvt WVTransform {mustBeNonempty}
                options.alpha (1,1) double {mustBeNonnegative} = 1/(200*86400) % linear bottom friction, try 1/(200*86400) https://www.nemo-ocean.eu/doc/node70.html
            end
            self@WVForcing(wvt,"thermal damping",WVForcingType("PVSpatial"));
            self.alpha = options.alpha;
            self.alpha_scaled = self.alpha/wvt.Lr2;
        end

        function Fpv = addPotentialVorticitySpatialForcing(self, wvt, Fpv)
            Fpv = Fpv + self.alpha_scaled * wvt.psi;
        end

        function force = forcingWithResolutionOfTransform(self,wvtX2)
            force = WVThermalDamping(wvtX2,alpha=self.alpha);
        end
    end
    methods (Static)
        function vars = classRequiredPropertyNames()
            vars = {'alpha'};
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            arguments (Output)
                propertyAnnotations CAPropertyAnnotation
            end
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
            propertyAnnotations(end+1) = CANumericProperty('alpha', {}, 's^{-1}','thermal damping coefficient');
        end
    end

end