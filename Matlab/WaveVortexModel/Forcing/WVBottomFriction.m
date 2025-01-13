classdef WVBottomFriction < WVForcing
    % Bottom friction
    %
    % Applies linear bottom friction to the flow, i.e., $$\frac{du}{dt} = -r*u$$.
    %
    % - Topic: Initializing
    % - Declaration: WVBottomFriction < [WVForcing](/classes/wvforcing/)
    properties
        wvt
        r
    end
    properties (Dependent)
        uv_damp
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
            self@WVForcing("bottom friction");
            self.wvt = wvt;
            self.doesHydrostaticSpatialForcing = true;
            self.doesNonhydrostaticSpatialForcing = true;
            self.r = options.r;
        end

        function [Fu, Fv, Feta] = addHydrostaticSpatialForcing(self, wvt, Fu, Fv, Feta)
            Fu(:,:,1) = Fu(:,:,1) - self.r*wvt.u(:,:,1);
            Fv(:,:,1) = Fv(:,:,1) - self.r*wvt.v(:,:,1);
        end

        function [Fu, Fv, Fw, Feta] = addNonhydrostaticSpatialForcing(self, wvt, Fu, Fv, Fw, Feta)
            Fu(:,:,1) = Fu(:,:,1) - self.r*wvt.u(:,:,1);
            Fv(:,:,1) = Fv(:,:,1) - self.r*wvt.v(:,:,1);
        end

        function writeToFile(self,group,wvt)
            arguments
                self WVForcing {mustBeNonempty}
                group NetCDFGroup {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end
            group.addAttribute('r',self.r)
        end

        function force = forcingWithResolutionOfTransform(self,wvtX2)
            force = WVBottomFriction(wvtX2,r=self.r);
        end

        function flag = isequal(self,other)
            arguments
                self WVForcing
                other WVForcing
            end
            flag = isequal@WVForcing(self,other);
            flag = flag & isequal(self.r, other.r);
        end
    end

    methods (Static)
        function force = forcingFromFile(group,wvt)
            arguments
                group NetCDFGroup {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end
            force = WVBottomFriction(wvt,r=group.attributes('r'));
        end
    end

end