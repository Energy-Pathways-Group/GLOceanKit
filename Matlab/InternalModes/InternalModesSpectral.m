classdef InternalModesSpectral < InternalModesBase
    properties (Access = public)
    end
    
    properties (Access = private)
        n                   % number of points used in the eigenvalue problem
        z_lobatto_grid      % z_in coordinated, on Chebyshev extrema/Lobatto grid
        rho_lobatto         % rho on the above z_lobatto_grid
        rho_cheb            % rho in spectral space
        N2_cheb             % N2 in spectral space
        N2_z_lobatto        % N2 on the above z_lobatto_grid
        Diff1               % single derivative in spectral space
        T, T_z, T_zz        % Chebyshev polys (and derivs) on the z_lobatto_grid
        
        doesOutputGridSpanDomain = 0    % if not, the output grid can't be used for normalization
        canFastTransformOutput = 0      % if yes (and it spans the domain), a fast transform can be used
        T_out                           % *function* handle that transforms from spectral to z_out
        z_out                           % output grid
    end
    
    properties (Dependent)
        rho_z
        rho_zz
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = InternalModesSpectral(rho, z_in, z_out, latitude)
            self@InternalModesBase(rho,z_in,z_out,latitude);
            
            self.n = length(self.z_lobatto_grid);
            self.rho_cheb = fct(self.rho_lobatto);
            self.Diff1 = (2/self.Lz)*ChebyshevDifferentiationMatrix( self.n );
            self.N2_cheb= -(self.g/self.rho0)*self.Diff1*self.rho_cheb;
            self.N2_z_lobatto = ifct(self.N2_cheb);
            [self.T,self.T_z,self.T_zz] = ChebyshevPolynomialsOnGrid( self.z_lobatto_grid, self.n );
            
            self.InitializeOutputTransformation(z_out);
            self.N2 = self.T_out(self.N2_cheb);
        end
        
        % Superclass calls this method upon initialization when it
        % determines that the input is given in gridded form. The goal is
        % to initialize z_lobatto_grid and rho_lobatto.
        function self = InitializeWithGrid(self, rho, z_in)
            if (z_in(2) - z_in(1)) > 0
                z_in = flip(z_in);
                rho = flip(rho);
            end
            
            z_in_norm_grid = ChebyshevPolynomialsOnGrid( z_in ); % z_in, normalized to [-1, 1] (-1 bottom, 1 top)
            if IsChebyshevGrid(z_in) == 1
                self.z_lobatto_grid = z_in;
                self.rho_lobatto = rho;
            else
                N_points = FindSmallestGridWithNoGaps(z_in_norm_grid);
                self.z_lobatto_grid = (self.Lz/2)*(cos(((0:N_points-1)')*pi/(N_points-1)) + 1) + min(z_in); % z, on a chebyshev grid
                self.rho_lobatto = interp1(z_in, rho, self.z_lobatto_grid, 'spline'); % rho, interpolated to that grid
            end
        end
        
        % Superclass calls this method upon initialization when it
        % determines that the input is given in functional form. The goal
        % is to initialize z_lobatto_grid and rho_lobatto.
        function self = InitializeWithFunction(self, rho, z_min, z_max, z_out)
            % We have an analytical function describing density, so lets
            % place points on a Chebyshev grid such that z_out is fully
            % resolved.
            if (z_out(2) - z_out(1)) > 0 % make z_out decreasing
                z_out = flip(z_out);
            end
            
            % Note that z_out might not span the whole domain, so we may
            % need to add z_min and z_max to the end points.
            if z_max > z_out(1)
                z_in = cat(1,z_max,z_out);
            end
            if z_min < z_out(end)
                z_in = cat(1,z_out,z_min);
            end
            
            if IsChebyshevGrid(z_in) == 1
                self.z_lobatto_grid = z_in;
            else
                z_in_norm_grid = ChebyshevPolynomialsOnGrid( z_in ); % z_in, normalized to [-1, 1] (-1 bottom, 1 top)
                N_points = FindSmallestGridWithNoGaps(z_in_norm_grid);
                self.z_lobatto_grid = (self.Lz/2)*(cos(((0:N_points-1)')*pi/(N_points-1)) + 1) + z_min; % z, on a chebyshev grid
            end
            
            self.rho_lobatto = rho(self.z_lobatto_grid);
        end
        
        % After the input variables have been initialized, this is used to
        % initialize the output transformation, T_out(f).
        function self = InitializeOutputTransformation(self, z_out)
            if (min(z_out) == min(self.z_lobatto_grid) && max(z_out) == max(self.z_lobatto_grid))
                self.doesOutputGridSpanDomain = 1;
            else
                self.doesOutputGridSpanDomain = 0;
            end
            
            if self.doesOutputGridSpanDomain == 1 && IsChebyshevGrid(z_out) == 1
                self.canFastTransformOutput = 1; % May require truncation or padding, however.
                if length(z_out) == self.n
                    self.T_out = @(f_cheb) ifct(f_cheb);
                elseif length(z_out) > self.n
                    self.T_out = @(f_cheb) ifct(cat(1,f_cheb,zeros(length(z_out)-self.n,1)));
                elseif length(z_out) < self.n
                    self.T_out = @(f_cheb) ifct(f_cheb(1:length(z_out)));
                end
            else
                x = (2/self.Lz)*(z_out-min(self.z_lobatto_grid)) - 1;
                t = acos(x);
                TT = zeros(length(z_out),self.n);
                for iPoly=0:(self.n-1)
                    TT(:,iPoly+1) = cos(iPoly*t);
                end
                self.T_out = @(f_cheb) TT*f_cheb;
            end
            
            self.z_out = z_out;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computed (dependent) properties
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function value = get.rho_z(self)
            value = self.T_out(self.Diff1 * self.rho_lobatto);
        end
        
        function value = get.rho_zz(self)
            value = self.T_out(self.Diff1 * self.Diff1 * self.rho_lobatto);
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
            A = (self.T_zz - k*k*eye(self.n)*self.T);
            B = diag((self.f0*self.f0-self.N2_z_lobatto)/self.g)*self.T;
            
            % Lower boundary is rigid, G=0
            A(self.n,:) = self.T(self.n,:);
            B(self.n,:) = 0;
            
            % G=0 or G_z = \frac{1}{h_j} G at the surface, depending on the BC
            if strcmp(self.upper_boundary, 'free_surface')
                % G_z = \frac{1}{h_j} G at the surface
                A(1,:) = self.T_z(1,:);
                B(1,:) = self.T(1,:);
            elseif strcmp(self.upper_boundary, 'rigid_lid')
                A(1,:) = self.T(1,:);
                B(1,:) = 0;
            end
            
            h_func = @(lambda) 1.0 ./ lambda;
            [F,G,h] = ModesFromGEP(self,A,B,h_func);
        end
        
        % Take matrices A and B from the generalized eigenvalue problem
        % (GEP) and returns F,G,h. The h_func parameter is a function that
        % returns the eigendepth, h, given eigenvalue lambda from the GEP.
        function [F,G,h] = ModesFromGEP(self,A,B,h_func)
            [V,D] = eig( A, B );
            
            [lambda, permutation] = sort(real(diag(D)),1,'ascend');
            G_cheb=V(:,permutation);
            h = h_func(lambda.');
            
            F = zeros(length(self.z_out),self.n);
            G = zeros(length(self.z_out),self.n);
            
            if self.doesOutputGridSpanDomain == 1
                for j=1:self.n
                    G(:,j) = self.T_out(G_cheb(:,j));
                    F(:,j) = h(j) * self.T_out(self.Diff1 * G_cheb(:,j));
                end
                [F,G] = self.NormalizeModes(F,G,self.z_out);
            else
                error('This normalization condition is not yet implemented!')
            end
        end
        
    end
    
end

function bool = IsChebyshevGrid(z_in)
% make sure the grid is monotonically decreasing
if (z_in(2) - z_in(1)) > 0
    z_in = flip(z_in);
end

z_norm = ChebyshevPolynomialsOnGrid( z_in );
N_points = length(z_in);
xi=(0:N_points-1)';
z_lobatto_grid = cos(xi*pi/(N_points-1));
z_diff = z_norm-z_lobatto_grid;
if max(abs(z_diff)) < 1e-6
    bool = 1;
else
    bool = 0;
end
end

% Want to create a chebyshev grid that never has two or more point between
% its points. If that makes sense.
function N_points = FindSmallestGridWithNoGaps(z_in)
n = ceil(log2(length(z_in)));
N_points = 2^n;
z_lobatto_grid_grid = cos(((0:N_points-1)')*pi/(N_points-1));

while( length(unique(interp1(z_lobatto_grid_grid,z_lobatto_grid_grid,z_in,'previous'))) ~= length(z_in) )
    n = n + 1;
    N_points = 2^n;
    z_lobatto_grid_grid = cos(((0:N_points-1)')*pi/(N_points-1));
end

end
