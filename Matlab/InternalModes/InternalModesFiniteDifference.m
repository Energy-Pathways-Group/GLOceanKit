classdef InternalModesFiniteDifference < InternalModesBase
    properties (Access = public)
        orderOfAccuracy = 4
        Nz
    end
    
    properties (Access = private)
        Diff1
        Diff2
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = InternalModesFiniteDifference(rho, z_in, z_out, latitude, orderOfAccuracy)
            if nargin < 5
                orderOfAccuracy = 4;
            end
            
            obj@InternalModesBase(rho,z_in,z_out,latitude);
            
            if orderOfAccuracy < 2
                obj.orderOfAccuracy = 2;
            elseif orderOfAccuracy > obj.Nz
                obj.orderOfAccuracy = obj.Nz;
            end
        end
        
        function obj = InitWithOrder( orderOfAccuracy )

            
            % The operator D is the first order, central difference operator
            obj.Diff1 = FiniteDifferenceMatrix(1, obj.z, 1, 1, orderOfAccuracy);
            
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
        
        function [F,G,h] = ModesAtWavenumber(obj, k )
            % The eigenvalue equation is,
            % G_{zz} - K^2 G = \frac{f_0^2 -N^2}{gh_j}G
            % A = \frac{g}{f_0^2 -N^2} \left( \partial_{zz} - K^2*I \right)
            % B = I
            A = diag(-obj.g./(obj.N2 - obj.f0*obj.f0)) * (obj.Diff2 - k*k*eye(obj.Nz));
            B = eye(obj.Nz);
            
            % G=0 at the bottom (we chose this BC when creating Diff2)
            A(1,:) = obj.Diff2(1,:);
            B(1,:) = 0;
            
            % G=0 or G_z = \frac{1}{h_j} G at the surface, depending on the BC
            if strcmp(obj.upper_boundary, 'free_surface')
                % G_z = \frac{1}{h_j} G at the surface
                A(end,:) = obj.Diff2(end,:);
            elseif strcmp(obj.upper_boundary, 'rigid_lid')
                % G=0 at the surface (note we chose this BC when creating Diff2)
                A(end,:) = obj.Diff2(end,:);
                B(end,end)=0;
            end
            
            [V,D] = eig( A, B );
            
            [lambda, permutation] = sort(real(diag(D)),1,'ascend');
            G=V(:,permutation);
            h = (1.0 ./ lambda).';
            F = h * obj.Diff1 * G;
            
            [F,G] = obj.NormalizeModes(F,G,obj.z);
        end
        
        function [F,G,h] = ModesAtFrequency(obj, omega )
            % The eigenvalue equation is,
            % G_{zz} - K^2 G = \frac{f_0^2 -N^2}{gh_j}G
            % A = \frac{g}{f_0^2 -N^2} \left( \partial_{zz} - K^2*I \right)
            % B = I
            A = diag( (obj.f0*obj.f0 - omega*omega)./(obj.N2 - omega*omega) ) * obj.Diff2;
            B = eye(obj.Nz);
            
            % G=0 at the bottom (we chose this BC when creating Diff2)
            A(1,:) = obj.Diff2(1,:);
            B(1,:) = 0;
            
            % G=0 or G_z = \frac{1}{h_j} G at the surface, depending on the BC
            if strcmp(obj.upper_boundary, 'free_surface')
            prefactor = (omega*omega - obj.f0*obj.f0)/obj.g;
            A(end,:) = prefactor*obj.Diff2(end,:);
            elseif strcmp(obj.upper_boundary, 'rigid_lid')
                % G=0 at the surface (note we chose this BC when creating Diff2)
                A(end,:) = obj.Diff2(end,:);
                B(end,end)=0;
            end
            
            [V,D] = eig( A, B );
            
            [lambda, permutation] = sort(real(diag(D)),1,'ascend');
            G=V(:,permutation);
            h = ((omega*omega - obj.f0*obj.f0)./(obj.g*lambda)).';
            F = h * obj.Diff1 * G;
            
            [F,G] = obj.NormalizeModes(F,G,obj.z);
        end
    end
    
end