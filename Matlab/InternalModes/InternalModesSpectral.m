classdef InternalModesSpectral < InternalModes
    properties (Access = public)
        orderOfAccuracy = 4
        Nz
    end
    
    properties (Access = private)
        Diff1
        Diff2
        rho_z_norm  % rho on a chebyshev grid
        rho_cheb
        N2_cheb
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = InternalModesSpectral(rho, z_in, z_out, latitude)
            self@InternalModes(rho,z_in,z_out,latitude);
            
            self.rho_cheb = fct(rho_z_norm);
            self.Diff1 = (2/self.L_z)*ChebyshevDifferentiationMatrix( N_points );
            self.N2_cheb= -(self.g/self.rho0)*self.Diff1*self.rho_cheb;
        end
        
        % Superclass calls this method upon initialization
        function self = InitializeWithGrid(self, rho, z_in)
            if (z_in(2) - z_in(1)) > 0
                z_in = flip(z_in);
                rho = flip(rho);
            end
            
            z_in_norm_grid = ChebyshevPolynomialsOnGrid( z_in ); % z_in, normalized to [-1, 1] (-1 bottom, 1 top)
            if IsChebyshevGrid(z_in) == 1
                self.rho_z_norm = rho;
            else
                N_points = FindSmallestGridWithNoGaps(z_in_norm_grid);
                z_cheb_grid = cos(((0:N_points-1)')*pi/(N_points-1)); % z, on a chebyshev grid
                self.rho_z_norm = interp1(z_in_norm_grid, rho, z_cheb_grid, 'spline'); % rho, interpolated to that grid
            end      
        end
        
        function self = InitialOutputTransformation(self, z_out)
            
        end
        
        function [F,G,h] = ModesAtWavenumber(self, k )
            % The eigenvalue equation is,
            % G_{zz} - K^2 G = \frac{f_0^2 -N^2}{gh_j}G
            % A = \frac{g}{f_0^2 -N^2} \left( \partial_{zz} - K^2*I \right)
            % B = I
            A = (self.T_zz - k*k*eye(self.N_ev)*self.T);
            B = diag((self.f0*self.f0-self.N2)/self.g)*self.T;
            
            % Lower boundary is rigid, G=0
            A(self.N_ev,:) = self.T(self.N_ev,:);
            B(self.N_ev,:) = 0;
            
            % G=0 or G_z = \frac{1}{h_j} G at the surface, depending on the BC
            if strcmp(self.upper_boundary, 'free_surface')
                % G_z = \frac{1}{h_j} G at the surface
                A(1,:) = self.T_z(1,:);
                B(1,:) = self.T(1,:);
            elseif strcmp(self.upper_boundary, 'rigid_lid')
                A(1,:) = self.T(1,:);
                B(1,:) = 0;
            end
            
            [V,D] = eig( A, B );
            
            [lambda, permutation] = sort(real(diag(D)),1,'ascend');
            G_cheb=V(:,permutation);
            h = (1.0 ./ lambda).';
            
            F = h * self.Diff1 * G;
            
            [F,G] = self.NormalizeModes(F,G,self.z);
        end
        
        function [F,G] = ProjectOntoOutputGrid(self,G_cheb,h)
            
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
z_cheb = cos(xi*pi/(N_points-1));
z_diff = z_norm-z_cheb;
if max(abs(z_diff)) < 1e-6
    bool = 1;
else
    bool = 0;
end
end

% Want to create a chebyshev grid that never has two or more point between
% its points. If that makes sense.
function n = FindSmallestGridWithNoGaps(z_in)
n = ceil(log2(length(z_in)));
N_points = 2^n;
z_cheb_grid = cos(((0:N_points-1)')*pi/(N_points-1));

while( length(unique(interp1(z_cheb_grid,z_cheb_grid,z_in,'previous'))) ~= length(z_in) )
    n = n + 1;
    N_points = 2^n;
    z_cheb_grid = cos(((0:N_points-1)')*pi/(N_points-1));
end

end
