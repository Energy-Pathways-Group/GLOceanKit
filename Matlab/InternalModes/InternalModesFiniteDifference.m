classdef InternalModesFiniteDifference < InternalModesBase
    properties (Access = public)
        rho
        N2
        Nz
    end
    
    properties (Dependent)
        rho_z
        rho_zz
    end
    
    properties (Access = private)
        doesOutputGridSpanDomain = 0    % if not, the output grid can't be used for normalization
        n                   % length of z_diff
        z_diff              % the z-grid used for differentiation
        rho_z_diff          % rho on the z_diff grid
        N2_z_diff           % N2 on the z_diff grid
        Diff1               % 1st derivative matrix, w/ 1st derivative boundaries
        Diff2               % 2nd derivative matrix, w/ BCs set by upperBoundary property
        
        T_out               % *function* handle that transforms from z_diff functions to z_out functions
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = InternalModesFiniteDifference(rho, z_in, z_out, latitude, orderOfAccuracy)
            if nargin < 5
                orderOfAccuracy = 4;
            end
            
            if orderOfAccuracy < 2
                orderOfAccuracy = 2;
            elseif orderOfAccuracy > length(z_in)
                orderOfAccuracy = length(z_in);
            end
            
            self@InternalModesBase(rho,z_in,z_out,latitude,orderOfAccuracy);
            
            self.n = length(self.z_diff);
            self.Diff1 = FiniteDifferenceMatrix(1, self.z_diff, 1, 1, self.orderOfAccuracy);
            self.N2_z_diff = -(self.g/self.rho0) * self.Diff1 * self.rho_z_diff;
            self.upperBoundaryDidChange(); % this sets Diff2          
            
            self.InitializeOutputTransformation(z_out);
            self.rho = self.T_out(self.rho_z_diff);
            self.N2 = self.T_out(self.N2_z_diff);
        end
        
        % Superclass calls this method upon initialization when it
        % determines that the input is given in gridded form. The goal is
        % to initialize z_diff and rho_z_diff.
        function self = InitializeWithGrid(self, rho, z_in)
            self.z_diff = z_in;
            self.rho_z_diff = rho;
        end
        
        % Superclass calls this method upon initialization when it
        % determines that the input is given in functional form. The goal
        % is to initialize z_diff and rho_z_diff.
        function self = InitializeWithFunction(self, rho, z_min, z_max, z_out)
            if (min(z_out) == z_min && max(z_out) == z_max)
                self.z_diff = z_out;
                self.rho_z_diff = rho(self.z_diff);
            else
                error('Other cases not yet implemented');
                % Eventually we may want to use stretched coordinates as a
                % default
            end
        end
        
        % After the input variables have been initialized, this is used to
        % initialize the output transformation, T_out(f).
        function self = InitializeOutputTransformation(self, z_out)
            if (min(z_out) == min(self.z_diff) && max(z_out) == max(self.z_diff))
                self.doesOutputGridSpanDomain = 1;
            else
                self.doesOutputGridSpanDomain = 0;
            end
            
            if isequal(self.z_diff,z_out)
                self.T_out = @(f_in) f_in;
            else
                error('Other cases not yet implemented');
                % want to interpolate onto the output grid
            end
        end
        
        % Superclass property setter is set to call this method after a
        % change.
        function self = upperBoundaryDidChange(self)
            if strcmp(self.upperBoundary, 'free_surface')
                rightBCDerivs = 1;
            else
                rightBCDerivs = 0;
            end
            self.Diff2 = FiniteDifferenceMatrix(2, self.z_diff, 0, rightBCDerivs, self.orderOfAccuracy);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computed (dependent) properties
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function value = get.rho_z(self)
            value = self.Diff1 * self.rho_z_diff;
        end
        
        function value = get.rho_zz(self)
            diff2 = FiniteDifferenceMatrix(2, self.z_diff, 2, 2, self.orderOfAccuracy);
            value = diff2 * self.rho_z_diff;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computation of the modes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [F,G,h] = ModesAtWavenumber(self, k )
            % The eigenvalue equation is,
            % G_{zz} - K^2 G = \frac{f_0^2 -N^2}{gh_j}G
            % A = \frac{g}{f_0^2 -N^2} \left( \partial_{zz} - K^2*I \right)
            % B = I
            A = diag(-self.g./(self.N2 - self.f0*self.f0)) * (self.Diff2 - k*k*eye(self.n));
            B = eye(self.n);
            
            % Bottom boundary condition (always taken to be G=0)
            % NOTE: we already chose the correct BCs when creating the
            % Diff2 matrix
            A(1,:) = self.Diff2(1,:);
            B(1,:) = 0;
            
            % Surface boundary condition
            A(end,:) = self.Diff2(end,:);
            if strcmp(self.upperBoundary, 'free_surface')
                % G_z = \frac{1}{h_j} G at the surface
                B(end,end)=1;
            elseif strcmp(self.upperBoundary, 'rigid_lid')
                % G=0 at the surface (note we chose this BC when creating Diff2)
                B(end,end)=0;
            end
            
            h_func = @(lambda) 1.0 ./ lambda;
            [F,G,h] = ModesFromGEP(self,A,B,h_func);
        end
        
        function [F,G,h] = ModesAtFrequency(self, omega )
            % The eigenvalue equation is,
            % G_{zz} - K^2 G = \frac{f_0^2 -N^2}{gh_j}G
            % A = \frac{g}{f_0^2 -N^2} \left( \partial_{zz} - K^2*I \right)
            % B = I
            
            %%%%% Do we get better accuracy by not dividing by N2??? 
            A = diag( (self.f0*self.f0 - omega*omega)./(self.N2 - omega*omega) ) * self.Diff2;
            B = eye(self.n);
            
            % Bottom boundary condition (always taken to be G=0)
            % NOTE: we already chose the correct BCs when creating the
            % Diff2 matrix
            A(1,:) = self.Diff2(1,:);
            B(1,:) = 0;
            
            % Surface boundary condition
            if strcmp(self.upperBoundary, 'free_surface')
                prefactor = (omega*omega - self.f0*self.f0)/self.g;
                A(end,:) = prefactor*self.Diff2(end,:);
                B(end,end)=1;
            elseif strcmp(self.upperBoundary, 'rigid_lid')
                % G=0 at the surface (note we chose this BC when creating Diff2)
                A(end,:) = self.Diff2(end,:);
                B(end,end)=0;
            end
            
            h_func = @(lambda) ((omega*omega - self.f0*self.f0)./(self.g*lambda)).';
            [F,G,h] = ModesFromGEP(self,A,B,h_func);
        end
        
        % Take matrices A and B from the generalized eigenvalue problem
        % (GEP) and returns F,G,h. The h_func parameter is a function that
        % returns the eigendepth, h, given eigenvalue lambda from the GEP.
        function [F,G,h] = ModesFromGEP(self,A,B,h_func)
            [V,D] = eig( A, B );
            
            [lambda, permutation] = sort(real(diag(D)),1,'ascend');
            G = V(:,permutation);
            h = h_func(lambda.');
            
            F = zeros(length(self.z),self.n);
            for j=1:self.n
                F(:,j) = h(j) * self.Diff1 * G(:,j);
            end
            
            if self.doesOutputGridSpanDomain == 1
                [F,G] = self.NormalizeModes(F,G,self.z);
            else
                error('This normalization condition is not yet implemented!')
            end
        end

    end
    
end