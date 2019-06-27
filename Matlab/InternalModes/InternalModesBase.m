classdef (Abstract) InternalModesBase < handle
    % InternalModesBase is an abstract class from which four different
    % internal mode solution methods derive.
    %
    % The class InternalModes is a wrapper around these four classes.
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
    %       [F,G,h,omega] = modes.ModesAtWavenumber(0.01);
    %   or frequency,
    %       [F,G,h,k] = modes.ModesAtWavenumber(5*modes.f0);
    %
    %   Jeffrey J. Early
    %   jeffrey@jeffreyearly.com
    %
    %   March 14th, 2017        Version 1.0
    %
    %   See also INTERNALMODES, INTERNALMODESSPECTRAL,
    %   INTERNALMODESDENSITYSPECTRAL, INTERNALMODESWKBSPECTRAL, and
    %   INTERNALMODESFINITEDIFFERENCE.

    
    properties (Access = public)
        shouldShowDiagnostics = 0 % flag to show diagnostic information, default = 0
        
        latitude % Latitude for which the modes are being computed.
        f0 % Coriolis parameter at the above latitude.
        Lz % Depth of the ocean.
        rho0 % Density at the surface of the ocean.
        
        nModes % Used to limit the number of modes to be returned.
        
        z % Depth coordinate grid used for all output (same as zOut).
        zDomain % [zMin zMax]
        requiresMonotonicDensity = 0
        
        gridFrequency = [] % last requested frequency from the user---set to f0 if a wavenumber was last requested
        normalization = Normalization.kConstant % Normalization used for the modes. Either Normalization.(kConstant, omegaConstant, uMax, or wMax).
        upperBoundary = UpperBoundary.rigidLid  % Surface boundary condition. Either UpperBoundary.rigidLid (default) or UpperBoundary.freeSurface.
        lowerBoundary = LowerBoundary.freeSlip  % Lower boundary condition. Either LowerBoundary.freeSlip (default) or LowerBoundary.none.
    end
    
    properties (Dependent)
        zMin
        zMax
    end
    
    properties (Abstract)
        rho  % Density on the z grid.
        N2 % Buoyancy frequency on the z grid, $N^2 = -\frac{g}{\rho(0)} \frac{\partial \rho}{\partial z}$.
        rho_z % First derivative of density on the z grid.
        rho_zz % Second derivative of density on the z grid.
    end
    
    properties (Access = protected)
        g = 9.81 % 9.81 meters per second.
        omegaFromK % function handle to compute omega(h,k)
        kFromOmega % function handle to compute k(h,omega)
    end
    
    methods (Abstract)
        [F,G,h,omega] = ModesAtWavenumber(self, k ) % Return the normal modes and eigenvalue at a given wavenumber.
        [F,G,h,k] = ModesAtFrequency(self, omega ) % Return the normal modes and eigenvalue at a given frequency.
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
            if  (norm ~= Normalization.kConstant && norm ~= Normalization.omegaConstant && norm ~= Normalization.uMax && norm ~= Normalization.wMax)
                error('Invalid normalization! Valid options: Normalization.kConstant, Normalization.omegaConstant, Normalization.uMax, Normalization.wMax')
            else
                obj.normalization = norm;
            end
        end
        
        function set.upperBoundary(self,upperBoundary)
            self.upperBoundary = upperBoundary;
            self.upperBoundaryDidChange();
        end
        
        function value = get.zMin(self)
            value = self.zDomain(1);
        end
        
        function value = get.zMax(self)
            value = self.zDomain(2);
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
            
            self.zDomain = [min(z_in) max(z_in)];
            self.Lz = self.zDomain(2)-self.zDomain(1);
            self.latitude = latitude;
            self.f0 = 2*(7.2921e-5)*sin(latitude*pi/180);
            self.z = z_out; % Note that z might now be a col-vector, when user asked for a row-vector.
            
            % Set properties supplied as name,value pairs
            userSpecifiedRho0 = 0;
            nargs = length(varargin);
            if mod(nargs,2) ~= 0
                error('Arguments must be given as name/value pairs');
            end
            for k = 1:2:length(varargin)
                self.(varargin{k}) = varargin{k+1};
                if strcmp(varargin{k}, 'rho0')
                    userSpecifiedRho0 = 1;
                end
            end
             
            % Is density specified as a function handle or as a grid of
            % values?
            if isa(rho,'function_handle') == true
                if numel(z_in) ~= 2
                    error('When using a function handle, z_domain must be an array with two values: z_domain = [z_bottom z_surface];')
                end
                if userSpecifiedRho0 == 0
                    self.rho0 = rho(max(z_in));
                end
                if self.shouldShowDiagnostics == 1
                    fprintf('Initialized %s class with a function handle.\n', class(self));
                end
                self.InitializeWithFunction(rho, min(z_in), max(z_in));
            elseif isa(rho,'numeric') == true
                if numel(rho) ~= length(rho) || length(rho) ~= length(z_in)
                    error('rho must be 1 dimensional and z must have the same length');
                end
                if isrow(rho)
                    rho = rho.';
                end
                if userSpecifiedRho0 == 0
                    self.rho0 = min(rho);
                end
                if self.shouldShowDiagnostics == 1
                    fprintf('Initialized %s class with gridded data.\n', class(self));
                end
                [z_in,I] = sort(z_in,'ascend');
                rho = rho(I);
                self.InitializeWithGrid(rho,z_in);
            else
                error('rho must be a function handle or an array.');
            end   
            
            self.kFromOmega = @(h,omega) sqrt((omega^2 - self.f0^2)./(self.g * h));
            self.omegaFromK = @(h,k) sqrt( self.g * h * k^2 + self.f0^2 );
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
        function [F,G,F2,N2G2] = NormalizeModes(self,F,G,N2,z)
            % This method normalizes the modes F,G using trapezoidal
            % integration on the given z grid. At the moment, this is only
            % used by the finite differencing algorithm, as the spectral
            % methods can use a superior (more accurate) technique of
            % directly integrating the polynomials.
            if z(2)-z(1) > 0
                direction = 'last';
            else
                direction = 'first';
            end
            
            N2G2 = zeros(1,length(G(1,:)));
            F2 = zeros(1,length(G(1,:)));
            
            [maxIndexZ] = find(N2-self.gridFrequency*self.gridFrequency>0,1,direction);  
            for j=1:length(G(1,:))
                switch self.normalization
                    case Normalization.uMax
                        A = max( abs(F(:,j)) );
                        G(:,j) = G(:,j) / A;
                        F(:,j) = F(:,j) / A;
                    case Normalization.wMax
                        A = max( abs(G(:,j)) );
                        G(:,j) = G(:,j) / A;
                        F(:,j) = F(:,j) / A;
                    case Normalization.kConstant
                        if z(2)-z(1) > 0
                            G20 = G(end,j)^2;
                        else
                            G20 = G(1,j)^2;
                        end
                        A = abs(G20 + trapz( z, (1/self.g) * (N2 - self.f0*self.f0) .* G(:,j) .^ 2));
                        G(:,j) = G(:,j) / sqrt(A);
                        F(:,j) = F(:,j) / sqrt(A);
                    case Normalization.omegaConstant
                        A = abs(trapz( z, (1/abs(z(end)-z(1))) .* F(:,j) .^ 2));
                        G(:,j) = G(:,j) / sqrt(A);
                        F(:,j) = F(:,j) / sqrt(A);
                end
                
                if F(maxIndexZ,j)< 0
                    F(:,j) = -F(:,j);
                    G(:,j) = -G(:,j);
                end
                
                F2(j) = abs(trapz( z, F(:,j) .^ 2));
                N2G2(j) = abs(trapz(z, N2.* (G(:,j).^2)));
            end
        end
        
    end
end

