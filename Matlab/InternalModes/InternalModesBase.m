classdef (Abstract) InternalModesBase < handle
    % InternalModesBase is an abstract class from which four different
    % internal mode solution methods derive.
    %
    % The class INTERNALMODES is a wrapper around these four classes.
    %
    %   There are two primary methods of initializing this class: either
    %   you specify density with gridded data, or you specify density as an
    %   analytical function.
    %
    %   modes = InternalModes(rho,z,zOut,latitude);
    %   'rho' must be a vector of gridded density data with length matching
    %   'z'. zOut is the output grid upon which all returned function will
    %   be given.
    %
    %   modes = InternalModes(rho,zDomain,zOut,latitude);
    %   'rho' must be a function handle, e.g.
    %       rho = @(z) -(N0*N0*rho0/g)*z + rho0
    %   zDomain must be an array with two values: z_domain = [z_bottom
    %   z_surface];
    %
    %   Once initialized, you can request variations of the density, e.g.,
    %       N2 = modes.N2;
    %       rho_zz = modes.rho_zz;
    %   or you can request the internal modes at a given wavenumber,
    %       [F,G,h] = modes.ModesAtWavenumber(0.01);
    %   or frequency,
    %       [F,G,h] = modes.ModesAtWavenumber(5*modes.f0);
    properties (Access = public)
        latitude % Latitude for which the modes are being computed.
        f0 % Coriolis parameter at the above latitude.
        Lz % Depth of the ocean.
        rho0 % Density at the surface of the ocean.
        
        nModes % Used to limit the number of modes to be returned.
        
        z % Depth coordinate grid used for all output (same as zOut).
        
        normalization = 'const_G_norm' % Normalization used for the modes. Either 'const_G_norm' (default), 'const_F_norm', 'max_u' or 'max_w'.
        upperBoundary = 'rigid_lid' % Surface boundary condition. Either 'rigid_lid' (default) or 'free_surface'.
    end
    
    properties (Abstract)
        rho  % Density on the z grid.
        N2 % Buoyancy frequency on the z grid, $N^2 = -\frac{g}{\rho(0)} \frac{\partial \rho}{\partial z}$.
        rho_z % First derivative of density on the z grid.
        rho_zz % Second derivative of density on the z grid.
    end
    
    properties (Access = protected)
        g = 9.81 % 9.81 meters per second.
    end
    
    methods (Abstract)
        [F,G,h] = ModesAtWavenumber(self, k ) % Return the normal modes and eigenvalue at a given wavenumber.
        [F,G,h] = ModesAtFrequency(self, omega ) % Return the normal modes and eigenvalue at a given frequency.
    end
    
    methods (Abstract, Access = protected)
        self = InitializeWithGrid(self, rho, z_in) % Used internally by subclasses to intialize with a density grid.
        self = InitializeWithFunction(self, rho, z_min, z_max, z_out) % Used internally by subclasses to intialize with a density function.
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Set the normalization and upper boundary condition
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set.normalization(obj,norm)
            if  (~strcmp(norm, 'max_u') && ~strcmp(norm, 'max_w') && ~strcmp(norm, 'const_G_norm') && ~strcmp(norm, 'const_F_norm'))
                error('Invalid normalization!')
            else
                obj.normalization = norm;
            end
        end
        
        function set.upperBoundary(self,upperBoundary)
            if  (~strcmp(upperBoundary, 'free_surface') && ~strcmp(upperBoundary, 'rigid_lid') )
                error('Invalid upper boundary condition!')
            else
                self.upperBoundary = upperBoundary;
                self.upperBoundaryDidChange();
            end
        end
    end
    
    methods (Access = protected)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = InternalModesBase(rho, z_in, z_out, latitude, varargin)
            % Initialize with either a grid or analytical profile.
            
            % Make everything a column vector
            if isrow(z_in)
                z_in = z_in.';
            end
            if isrow(z_out)
                z_out = z_out.';
            end
            
            self.Lz = max(z_in) - min(z_in);
            self.latitude = latitude;
            self.f0 = 2*(7.2921e-5)*sin(latitude*pi/180);
            self.z = z_out; % Note that z might now be a col-vector, when user asked for a row-vector.
            
            % Set properties supplied as name,value pairs
            nargs = length(varargin);
            if mod(nargs,2) ~= 0
                error('Arguments must be given as name/value pairs');
            end
            for k = 1:2:length(varargin)
                self.(varargin{k}) = varargin{k+1};
            end
             
            % Is density specified as a function handle or as a grid of
            % values?
            if isa(rho,'function_handle') == true
                if numel(z_in) ~= 2
                    error('When using a function handle, z_domain must be an array with two values: z_domain = [z_bottom z_surface];')
                end
                self.rho0 = rho(max(z_in));
                self.InitializeWithFunction(rho, min(z_in), max(z_in), z_out);
            elseif isa(rho,'numeric') == true
                if numel(rho) ~= length(rho) || length(rho) ~= length(z_in)
                    error('rho must be 1 dimensional and z must have the same length');
                end
                if isrow(rho)
                    rho = rho.';
                end
                self.rho0 = min(rho);
                self.InitializeWithGrid(rho,z_in);
            else
                error('rho must be a function handle or an array.');
            end   
        end
        

        
        function self = upperBoundaryDidChange(self)
            % This function is called when the user changes the surface
            % boundary condition. By overriding this function, a subclass
            % can respond as necessary.
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Generical function to normalize
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [F,G] = NormalizeModes(self,F,G,N2,z)
            % This method normalizes the modes F,G using trapezoidal
            % integration on the given z grid. At the moment, this is only
            % used by the finite differencing algorithm, as the spectral
            % methods can use a superior (more accurate) technique of
            % directly integrating the polynomials.
            [~,maxIndexZ] = max(z);
            for j=1:length(G(1,:))
                if strcmp(self.normalization, 'max_u')
                    A = max( abs(F(:,j)) );
                    G(:,j) = G(:,j) / A;
                    F(:,j) = F(:,j) / A;
                elseif strcmp(self.normalization, 'max_w')
                    A = max( abs(G(:,j)) );
                    G(:,j) = G(:,j) / A;
                    F(:,j) = F(:,j) / A;
                elseif strcmp(self.normalization, 'const_G_norm')
                    A = abs(trapz( z, (1/self.g) * (N2 - self.f0*self.f0) .* G(:,j) .^ 2));
                    G(:,j) = G(:,j) / sqrt(A);
                    F(:,j) = F(:,j) / sqrt(A);
                elseif strcmp(self.normalization, 'const_F_norm')
                    A = abs(trapz( z, (1/abs(z(end)-z(1))) .* F(:,j) .^ 2));
                    G(:,j) = G(:,j) / sqrt(A);
                    F(:,j) = F(:,j) / sqrt(A);
                end
                
                if F(maxIndexZ,j)< 0
                    F(:,j) = -F(:,j);
                    G(:,j) = -G(:,j);
                end
            end
        end
        
    end
end

