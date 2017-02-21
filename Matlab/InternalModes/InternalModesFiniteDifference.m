classdef InternalModesFiniteDifference < InternalModes
    properties (Access = public)
        orderOfAccuracy = 4
    end
    
    properties (Access = private)
        Diff1
        Diff2
    end
    
    methods
        function obj = InitWithOrder( orderOfAccuracy )
            if orderOfAccuracy < 2
                obj.orderOfAccuracy = 2;
            elseif orderOfAccuracy > N
                obj.orderOfAccuracy = N;
            end
            
            % The operator D is the first order, central difference operator
            obj.Diff1 = FiniteDifferenceMatrix(1, obj.z, 1, 1, orderOfAccuracy);
            
            %rho_0=trapz( z, rho)/(z(end)-z(1));
            rho_0=min(rho);
            
            % Convert density to buoyancy.
            Bouyancy=-g*rho/rho_0;
            
            % N^2 is computed assuming that d^2rho/dz^2=0 at the end points.
            obj.N2=obj.Diff1*Bouyancy;
            
            if strcmp(obj.upper_boundary, 'free_surface')
                rightBCDerivs = 1;
            else
                rightBCDerivs = 0;
            end
            obj.Diff2 = FiniteDifferenceMatrix(2, obj.z, 0, rightBCDerivs, orderOfAccuracy);
        end
        
        function obj = ModesAtWavenumber(obj, k )
            
        end
    end
    
end