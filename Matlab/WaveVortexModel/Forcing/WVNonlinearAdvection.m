classdef WVNonlinearAdvection < WVForcing
    % The advective flux, $$u\cdot \nabla u$$
    %
    % This computes the nonlinear advection terms.
    %
    % - Topic: Initializing
    % - Declaration: WVAdvectiveFlux < [WVForcing](/classes/wvforcing/)
    properties
        dLnN2 = 0
    end

    methods
        function self = WVNonlinearAdvection(wvt)
            % initialize the WVAdvectiveFlux nonlinear flux
            %
            % - Declaration: self = WVAdvectiveFlux(wvt,options)
            % - Parameter wvt: a WVTransform instance
            arguments
                wvt WVTransform {mustBeNonempty}
            end
            self@WVForcing(wvt,"nonlinear advection",WVForcingType(["HydrostaticSpatial" "NonhydrostaticSpatial" "PVSpatial"]));
            self.priority = 127;
            if isa(wvt,'WVStratification')
                self.dLnN2 = shiftdim(wvt.dLnN2,-2);
            end
        end
        
        function [Fu, Fv, Feta] = addHydrostaticSpatialForcing(self, wvt, Fu, Fv, Feta)
            Fu = Fu - (wvt.u .* wvt.diffX(wvt.u)   + wvt.v .* wvt.diffY(wvt.u)   + wvt.w .*  wvt.diffZF(wvt.u));
            Fv = Fv - (wvt.u .* wvt.diffX(wvt.v)   + wvt.v .* wvt.diffY(wvt.v)   + wvt.w .*  wvt.diffZF(wvt.v));
            Feta = Feta - (wvt.u .* wvt.diffX(wvt.eta) + wvt.v .* wvt.diffY(wvt.eta) + wvt.w .* (wvt.diffZG(wvt.eta) + wvt.eta .* self.dLnN2));
        end

        function [Fu, Fv, Fw, Feta] = addNonhydrostaticSpatialForcing(self, wvt, Fu, Fv, Fw, Feta)
            Fu = Fu - (wvt.u .* wvt.diffX(wvt.u)   + wvt.v .* wvt.diffY(wvt.u)   + wvt.w .*  wvt.diffZF(wvt.u));
            Fv = Fv - (wvt.u .* wvt.diffX(wvt.v)   + wvt.v .* wvt.diffY(wvt.v)   + wvt.w .*  wvt.diffZF(wvt.v));
            Fw = Fw - (wvt.u .* wvt.diffX(wvt.w)   + wvt.v .* wvt.diffY(wvt.w)   + wvt.w .*  wvt.diffZG(wvt.w));
            Feta = Feta - (wvt.u .* wvt.diffX(wvt.eta) + wvt.v .* wvt.diffY(wvt.eta) + wvt.w .* (wvt.diffZG(wvt.eta) + wvt.eta .* self.dLnN2));
        end

        function Fpv = addPotentialVorticitySpatialForcing(self, wvt, Fpv)
            Fpv = Fpv - (wvt.u.*wvt.diffX(wvt.qgpv) + wvt.v.*wvt.diffY(wvt.qgpv));
        end

        function force = forcingWithResolutionOfTransform(self,wvtX2)
            force = WVNonlinearAdvection(wvtX2);
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