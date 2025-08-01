classdef WVBetaPlanePVAdvection < WVForcing
    % Advection of QGPV from beta
    %
    % This applies \beta*v_g to the PV (A0) flux of a simulation. Note that
    % this may not be justified for Boussinesq flow, but it works.
    %
    %
    % - Topic: Initializing
    % - Declaration: WVBetaPlanePVAdvection < [WVForcing](/classes/wvforcing/)
    properties
        betaA0
    end
    methods
        function self = WVBetaPlanePVAdvection(wvt)
            arguments
                wvt WVTransform {mustBeNonempty}
            end
            self@WVForcing(wvt,"beta-plane advection of qgpv",WVForcingType(["Spectral" "PVSpatial"]));
            self.betaA0 = -wvt.beta * wvt.VA0;
            self.betaA0(1,1,1) = 0;
        end

        function Fpv = addPotentialVorticitySpatialForcing(self, wvt, Fpv)
            Fpv = Fpv - wvt.beta * wvt.v;
        end

        function [Fp, Fm, F0] = addSpectralForcing(self, wvt, Fp, Fm, F0)
            F0 = F0 + self.betaA0 .* wvt.A0;
        end

        function force = forcingWithResolutionOfTransform(self,wvtX2)
            force = WVBetaPlanePVAdvection(wvtX2);
        end
    end

    methods (Static)
        function vars = classRequiredPropertyNames()
            vars = {};
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            arguments (Output)
                propertyAnnotations CAPropertyAnnotation
            end
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
        end
    end

end