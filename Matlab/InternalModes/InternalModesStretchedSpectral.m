classdef InternalModesStretchedSpectral < InternalModesSpectral
    properties %(Access = private)
        
        s_lobatto_grid      % stretched density coordinate, on Chebyshev extrema/Lobatto grid
        z_s_lobatto         % The value of z, at the s_lobatto_grid points
        T_s                 % *function* handle that transforms from stretched-spectral to z_s_lobatto
        
        s_out               % desired locations of the output in s-coordinate (deduced from z_out)
        
        N2_s_grid           % N2 on the z_s_lobatto grid
        N2_z_s_grid         % (d/dz)N2 on the z_s_lobatto grid
    end
    

    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = InternalModesStretchedSpectral(rho, z_in, z_out, latitude)
            self@InternalModesSpectral(rho,z_in,z_out,latitude);
            
            % Need InitializeChebyshevMatrices to be called with
            % s_lobatto_grid and s_out
            
            % The eigenvalue problem will be solved using N2 and N2_z, so
            % now we need transformations to project them onto the
            % stretched grid
            t = acos((2/self.Lz)*(self.z_s_lobatto-min(self.z_lobatto_grid)) - 1);
            self.T_s = zeros(length(t),self.n);
            for iPoly=0:(self.n-1)
                self.T_s(:,iPoly+1) = cos(iPoly*t);
            end
            
            self.N2_s_grid = self.T_s * (self.N2_cheb);
            self.N2_z_s_grid = self.T_s * (self.Diff1*self.N2_cheb);
        end
        
        % Superclass calls this method upon initialization when it
        % determines that the input is given in gridded form. The goal is
        % to initialize z_lobatto_grid and rho_lobatto.
        function self = InitializeWithGrid(self, rho, z_in)
            % Superclass initializes z_lobatto_grid and rho_lobatto
            InitializeWithGrid@InternalModesSpectral(self, rho, z_in);
            
            s = -self.g*rho/self.rho0;
            L_s = max(s)-min(s);
            N_points = length(z_in); % More is better, but this should be sufficient.
            xi=(0:N_points-1)';
            self.s_lobatto_grid = abs(L_s/2)*(cos(xi*pi/(N_points-1))+1) + min(s); % s Chebyshev grid on the s-coordinate
            self.z_s_lobatto = interp1(s, z_in, self.s_lobatto_grid, 'spline'); % z, evaluated on that s grid
            
            self.s_out = interp1(z_in, s, self.z, 'spline'); % z, evaluated on that s grid
        end
        
        % Superclass calls this method upon initialization when it
        % determines that the input is given in functional form. The goal
        % is to initialize z_lobatto_grid and rho_lobatto.
        function self = InitializeWithFunction(self, rho, z_min, z_max, z_out)
            % Superclass initializes z_lobatto_grid and rho_lobatto
            InitializeWithFunction@InternalModesSpectral(self, rho, z_min, z_max, z_out);
            
            % Create a stretched grid that includes all the points of z_out
            s = @(z) -self.g*rho(z)/self.rho0;
            s_min = s(z_max);
            s_max = s(z_min);
            Ls = s_max-s_min;
            
            s_grid = s(z_out);
            if (s_grid(2) - s_grid(1)) > 0 % make s_grid decreasing
                s_grid = flip(s_grid);
            end
            if s_max > s_grid(1)
                s_grid = cat(1,s_max,s_grid);
            end
            if s_min < s_grid(end)
                s_grid = cat(1,s_grid,s_min);
            end
            
            if self.IsChebyshevGrid(s_grid) == 1
                self.s_lobatto_grid = s_grid;
            else
                % This will likely be waaay overkill in most cases... need
                % smarter logic here.
                N_points = self.FindSmallestGridWithNoGaps(s_grid);
                self.s_lobatto_grid = (Ls/2)*(cos(((0:N_points-1)')*pi/(N_points-1)) + 1) + s_min; % z, on a chebyshev grid
            end
            
            % Now create a transformation for functions defined on
            % z_lobatto to (spectrally) take them into s_lobatto.
            % We use the fact that we have a function handle to iteratively
            % improve this projection.
            self.z_s_lobatto = interp1(s(self.z_lobatto_grid), self.z_lobatto_grid, self.s_lobatto_grid, 'spline');
            for i=1:5
                self.z_s_lobatto = interp1(s(self.z_s_lobatto), self.z_s_lobatto, self.s_lobatto_grid, 'spline');
                fprintf('max diff: %f\n', max(s(self.z_s_lobatto) - self.s_lobatto_grid)/max(self.s_lobatto_grid));
            end
            
            self.s_out = s(z_out);
        end
        
        function self = InitializeChebyshevTMatrices(self)
            self.InitializeChebyshevTMatricesCallout(self.s_lobatto_grid,self.s_out);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computation of the modes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [F,G,h] = ModesAtWavenumber(self, k )
            A = diag(self.N2_s_grid .* self.N2_s_grid)*self.T_zz + diag(self.N2_z_s_grid)*self.T_z - k*k*self.T;
            B = diag( (self.f0*self.f0 - self.N2_s_grid)/self.g )*self.T;
            
            % Lower boundary is rigid, G=0
            A(self.n,:) = self.T(self.n,:);
            B(self.n,:) = 0;
            
            % G=0 or G_z = \frac{1}{h_j} G at the surface, depending on the BC
            if strcmp(self.upperBoundary, 'free_surface')
                % G_z = \frac{1}{h_j} G at the surface
                A(1,:) = self.T_z(1,:);
                B(1,:) = self.T(1,:);
            elseif strcmp(self.upperBoundary, 'rigid_lid')
                A(1,:) = self.T(1,:);
                B(1,:) = 0;
            end
            
            h_func = @(lambda) 1.0 ./ lambda;
            [F,G,h] = self.ModesFromGEP(A,B,h_func);
        end
        
    end
    
end
