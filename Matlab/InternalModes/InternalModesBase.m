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
        
        rotationRate % rotation rate of the planetary body
        latitude % Latitude for which the modes are being computed.
        f0 % Coriolis parameter at the above latitude.
        Lz % Depth of the ocean.
        rho0 % Density at the surface of the ocean.
        
        nModes % Used to limit the number of modes to be returned.
        
        z % Depth coordinate grid used for all output (same as zOut).
        zDomain % [zMin zMax]
        requiresMonotonicDensity
        
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
        g % 9.81 meters per second.
        omegaFromK % function handle to compute omega(h,k)
        kFromOmega % function handle to compute k(h,omega)
    end
    
    methods (Abstract)
        [F,G,h,omega] = ModesAtWavenumber(self, k ) % Return the normal modes and eigenvalue at a given wavenumber.
        [F,G,h,k] = ModesAtFrequency(self, omega ) % Return the normal modes and eigenvalue at a given frequency.
    end
    
    methods (Abstract, Access = protected)
        self = InitializeWithBSpline(self, rho) % Used internally by subclasses to intialize with a bspline.
        self = InitializeWithGrid(self, rho, z_in) % Used internally by subclasses to intialize with a density grid.
        self = InitializeWithFunction(self, rho, z_min, z_max, z_out) % Used internally by subclasses to intialize with a density function.
        self = InitializeWithN2Function(self, N2, z_min, z_max, z_out) % Used internally by subclasses to intialize with a density function.
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Set the normalization and upper boundary condition
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set.normalization(obj,norm)
            if  (norm ~= Normalization.geostrophic && norm ~= Normalization.kConstant && norm ~= Normalization.omegaConstant && norm ~= Normalization.uMax && norm ~= Normalization.wMax && norm ~= Normalization.surfacePressure)
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
        function self = InternalModesBase(options)
            arguments
                options.rho = ''
                options.N2 function_handle = @disp
                options.zIn (:,1) double = []
                options.zOut (:,1) double = []
                options.latitude (1,1) double = 33
                options.rho0 (1,1) double {mustBePositive} = 1025
                options.nModes (1,1) double = 0
                options.rotationRate (1,1) double = 7.2921e-5
                options.g (1,1) double = 9.81
            end
            if isempty(options.zIn)
                error('You must specify zIn');
            end

            self.requiresMonotonicDensity = self.requiresMonotonicDensitySetting();
            self.zDomain = [min(options.zIn) max(options.zIn)];
            self.Lz = self.zDomain(2)-self.zDomain(1);
            self.latitude = options.latitude;
            self.rotationRate = options.rotationRate;
            self.f0 = 2*(self.rotationRate)*sin(self.latitude*pi/180);
            self.rho0 = options.rho0;
            self.nModes = options.nModes;
            self.g = options.g;

            if isempty(options.zOut)
                self.z = options.zIn;
            else
                self.z = options.zOut;
            end
            
            if ~isequal(options.N2,@disp)
                self.InitializeWithN2Function(options.N2, min(options.zIn), max(options.zIn));
            elseif isa(options.rho,'function_handle') == true
                self.InitializeWithFunction(options.rho, min(options.zIn), max(options.zIn));
            elseif isa(options.rho,'BSpline') == true
                self.rho0 = rho(max(options.rho.domain));
                self.InitializeWithBSpline(options.rho);
            elseif isa(options.rho,'numeric') == true
                if length(options.rho) ~= length(options.zIn)
                    error('rho must be 1 dimensional and z must have the same length');
                end
                self.rho0 = min(options.rho);
                options.rho = reshape(options.rho,[],1);
                [zGrid,I] = sort(options.zIn,'ascend');
                rhoGrid = options.rho(I);
                self.InitializeWithGrid(rhoGrid,zGrid);
            else
                error('You must initialize InternalModes with rho, N2, rhoGrid, or rhoSpline');
            end
            
            self.kFromOmega = @(h,omega) sqrt((omega^2 - self.f0^2)./(self.g * h));
            self.omegaFromK = @(h,k) sqrt( self.g * h * k^2 + self.f0^2 );
        end
                
        function out=requiresMonotonicDensitySetting(~)
            out=0;
        end

        function self = upperBoundaryDidChange(self)
            % This function is called when the user changes the surface
            % boundary condition. By overriding this function, a subclass
            % can respond as necessary.
        end
   
    end
end

