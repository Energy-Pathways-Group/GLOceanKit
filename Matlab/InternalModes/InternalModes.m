
classdef InternalModes < handle
    % InternalModes This class solves the internal mode (Sturm-Liouville)
    % problem for a given density profile.
    %
    %   This computes the internal modes for u & v (F), the internal modes for w
    %   (G), the eigendepths (h). You can also request the buoyancy
    %   frequency (N2) and derivatives of rho (rho_z and rho_zz).
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
    %   or you can request the internal modes at a given wavenumber, k,
    %   where k is 2*pi/wavelength.
    %       [F,G,h] = modes.ModesAtWavenumber(0.01);
    %   or frequency,
    %       [F,G,h] = modes.ModesAtWavenumber(5*modes.f0);
    %
    %   There are two convenience methods,
    %       modes.ShowLowestModesAtWavenumber(0.0);
    %   and
    %       modes.ShowLowestModesAtFrequency(5*modes.f0)
    %   that can be used to quickly visual the modes.
    %
    %   Internally the InternalModes class is actually initializing one of
    %   four different classes used to solve the eigenvalue problem. By
    %   default is uses Cheyshev polynomials on a WKB coordinate grid. You
    %   can change this default by passing the name/value pair
    %   'method'/method where the method is either 'wkbSpectral',
    %   'densitySpectral', 'spectral' or 'finiteDifference'. For example,
    %       modes = InternalModes(rho,zDomain,zOut,lat, 'method',
    %       'finiteDifference');
    %   will initialize the class that uses finite differencing to solve
    %   the EVP.
    %
    %   You can also pass the argument 'nModes'/nModes to limit the number
    %   of modes that are returned.
    %
    %   Any other name value pairs are passed directly to the class being
    %   initialized.
    %
    %   This wrapper class also supports a number of basic test cases. You
    %   can initialize the class as follows,
    %       modes = InternalModes(stratification, method, n);
    %   where stratification must be either 'constant' (default) or
    %   'exponential', method must be one of the methods described above,
    %   and n is the number of points to be used.
    %
    %   When initialized with the built-in stratification profiles, you can
    %   use the functions,
    %       modes.ShowRelativeErrorAtWavenumber(0.01)
    %   or
    %       modes.ShowRelativeErrorAtFrequency(5*modes.f0)
    %   to assess the quality of the numerical methods. These test
    %   functions work for both constant G and F normalizations, as well as
    %   both rigid lid and free surface boundary conditions.
    %
    %   See also INTERNALMODESSPECTRAL, INTERNALMODESDENSITYSPECTRAL,
    %   INTERNALMODESWKBSPECTRAL, and INTERNALMODESFINITEDIFFERENCE
    %
    %
    %   Jeffrey J. Early
    %   jeffrey@jeffreyearly.com
    %
    %   March 14th, 2017        Version 1.0
    
    properties (Access = public)
        method % Numerical method used to solve the Sturm-Liouville equation. Either, 'wkbSpectral' (default), 'densitySpectral', 'spectral' or 'finiteDifference'.
        internalModes % Instance of actual internal modes class that is doing all the work.
    end
    
    properties (Dependent) 
        latitude % Latitude for which the modes are being computed.
        f0 % Coriolis parameter at the above latitude.
        nModes % Number of modes to be returned.
        
        Lz % Depth of the ocean.
        z % Depth coordinate grid used for all output (same as zOut).
        rho % Density on the z grid.
        N2 % Buoyancy frequency on the z grid, $N^2 = -\frac{g}{\rho(0)} \frac{\partial \rho}{\partial z}$.
        rho_z % First derivative of density on the z grid.
        rho_zz % Second derivative of density on the z grid.
        
        upperBoundary % Surface boundary condition. Either 'rigid_lid' (default) or 'free_surface'.
        normalization % Normalization used for the modes. Either 'const_G_norm' (default), 'const_F_norm', 'max_u' or 'max_w'.
    end
    
    properties (Access = private)
        isRunningTestCase = 0;
        stratification = 'user specified';
        rhoFunction
        N2Function
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function self = InternalModes(varargin)    
            % Initialize with either a grid or analytical profile.
            self.method = 'wkbSpectral';
            userSpecifiedMethod = 0;
            
            % First check to see if the user specified some extra arguments
            if nargin >= 5 
                extraargs = varargin(5:end);
                if mod(length(extraargs),2) ~= 0
                    error('Arguments must be given as name/value pairs.');
                end
                
                % now see if one of those arguments was specifying the
                % method
                for k = 1:2:length(extraargs)
                    if strcmp(extraargs{k}, 'method')
                        self.method = extraargs{k+1};
                        extraargs(k+1) = [];
                        extraargs(k) = [];
                        userSpecifiedMethod = 1;
                        break;
                    end
                end      
            else
                extraargs = {};
            end
            
            if nargin >= 4
                rho = varargin{1};
                zIn = varargin{2};
                zOut = varargin{3};
                latitude = varargin{4};
                
                isStratificationConstant = InternalModesConstantStratification.IsStratificationConstant(rho,zIn);
                
                if userSpecifiedMethod == 0 && isStratificationConstant == 1
                    fprintf('Initialization detected that you are using constant stratification. The modes will now be computed using the analytical form. If you would like to override this behavior, specify the method parameter.\n');
                    self.internalModes = InternalModesConstantStratification(rho,zIn,zOut,latitude,extraargs{:});
                elseif  strcmp(self.method, 'densitySpectral')
                    self.internalModes = InternalModesDensitySpectral(rho,zIn,zOut,latitude,extraargs{:});
                elseif  strcmp(self.method, 'wkbSpectral')
                    self.internalModes = InternalModesWKBSpectral(rho,zIn,zOut,latitude,extraargs{:});
                elseif strcmp(self.method, 'finiteDifference')
                    self.internalModes = InternalModesFiniteDifference(rho,zIn,zOut,latitude,extraargs{:});
                elseif strcmp(self.method, 'spectral')
                    self.internalModes = InternalModesSpectral(rho,zIn,zOut,latitude,extraargs{:});
                elseif strcmp(self.method, 'wkb')
                    self.internalModes = InternalModesWKB(rho,zIn,zOut,latitude,extraargs{:});
                else
                    error('Invalid method!')
                end
            else
                
                if nargin == 4
                   error('Invalid initialization'); 
                end
                if nargin < 3
                    n = 64;
                else
                    n = varargin{3};
                end
                if nargin < 2
                    theMethod = 'wkbSpectral';
                else
                    theMethod = varargin{2};
                end
                if nargin < 1
                    strat = 'constant';
                else
                    strat = varargin{1};
                end

                self.InitTestCase(strat,theMethod,n);
            end
            
        end
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Useful functions
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function ShowLowestModesAtWavenumber( self, k )
            % Quickly visualize the lowest modes at a given wavenumber.
            [F,G,h] = self.internalModes.ModesAtWavenumber( k );
            self.ShowLowestModesFigure(F,G,h);
        end
        
        function ShowLowestModesAtFrequency( self, omega )
            % Quickly visualize the lowest modes at a given frequency.
            [F,G,h] = self.internalModes.ModesAtFrequency( omega );
            self.ShowLowestModesFigure(F,G,h);
        end
        
        function ShowRelativeErrorAtWavenumber( self, k )
            % When using one of the built-in test cases, this method will
            % create a figure showing the relative error of the modes.
            if self.isRunningTestCase == 0
                error('Cannot show relative error for user specified stratification.\n');
            end
            
            [F,G,h] = self.internalModes.ModesAtWavenumber( k );
            
            % y is the true solution, x is the approximated
            errorFunction = @(x,y) max(abs(x-y),[],1)./max(abs(y),[],1);
            
            if  strcmp(self.stratification, 'constant')
                imConstant = InternalModesConstantStratification(self.rhoFunction,[-5000 0],self.z,self.latitude,'nModes',self.nModes);
                imConstant.upperBoundary = self.upperBoundary;
                imConstant.normalization = self.normalization;
                [F_analytical,G_analytical,h_analytical] = imConstant.ModesAtWavenumber( k );
            elseif strcmp(self.stratification, 'exponential')
                imAnalytical = InternalModesWKBSpectral(self.rhoFunction,[-5000 0],self.z,self.latitude,'nEVP',512,'nModes',self.nModes);
                imAnalytical.upperBoundary = self.upperBoundary;
                imAnalytical.normalization = self.normalization;
                [F_analytical,G_analytical,h_analytical] = imAnalytical.ModesAtWavenumber( k );
            else
                error('Invalid choice of stratification: you must use constant or exponential');
            end
            
            h_error = errorFunction(h,h_analytical);
            F_error = errorFunction(F,F_analytical);
            G_error = errorFunction(G,G_analytical);
            
            title = sprintf('Relative error of internal modes at k=%.2g with %d grid points using %s',k,length(self.z),self.fullMethodName);
            self.ShowErrorFigure(h_error,F_error,G_error,title);
        end
        
        function ShowRelativeErrorAtFrequency( self, omega )
            % When using one of the built-in test cases, this method will
            % create a figure showing the relative error of the modes.
            if self.isRunningTestCase == 0
                error('Cannot show relative error for user specified stratification.\n');
            end
            
            [F,G,h] = self.internalModes.ModesAtFrequency( omega );
            
            % y is the true solution, x is the approximated
            errorFunction = @(x,y) max(max(abs(x-y),[],1)./max(abs(y),[],1),1e-15);
            
            if  strcmp(self.stratification, 'constant')
                imConstant = InternalModesConstantStratification(self.rhoFunction,[-5000 0],self.z,self.latitude,'nModes',self.nModes);
                imConstant.upperBoundary = self.upperBoundary;
                imConstant.normalization = self.normalization;
                [F_analytical,G_analytical,h_analytical] = imConstant.ModesAtFrequency( omega );
            elseif strcmp(self.stratification, 'exponential')
                imAnalytical = InternalModesWKBSpectral(self.rhoFunction,[-5000 0],self.z,self.latitude,'nEVP',512,'nModes',self.nModes);
                imAnalytical.upperBoundary = self.upperBoundary;
                imAnalytical.normalization = self.normalization;
                [F_analytical,G_analytical,h_analytical] = imAnalytical.ModesAtFrequency( omega );
            else
                error('Invalid choice of stratification: you must use constant or exponential');
            end
            
            h_error = errorFunction(h,h_analytical);
            F_error = errorFunction(F,F_analytical);
            G_error = errorFunction(G,G_analytical);
            
            title = sprintf('Relative error of internal modes at omega=%.2gf0 with %d grid points using %s',omega/self.f0,length(self.z),self.fullMethodName);
            self.ShowErrorFigure(h_error,F_error,G_error,title);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Useful problem constants, latitude and f0
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set.latitude(self,value)
            self.internalModes.latitude = value;
        end
        function value = get.latitude(self)
            value = self.internalModes.latitude;
        end
        
        function set.f0(self,value)
            self.internalModes.f0 = value;
        end
        function value = get.f0(self)
            value = self.internalModes.f0;
        end
        
        function set.nModes(self,value)
            self.internalModes.nModes = value;
        end
        function value = get.nModes(self)
            value = self.internalModes.nModes;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Vertical grid and derivatives of the density profile on that grid
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function set.Lz(~,~)
            error('This property is readonly.')
        end
        function value = get.Lz(self)
            value = self.internalModes.Lz;
        end
        
        function set.z(~,~)
            error('This property is readonly.')
        end
        function value = get.z(self)
            value = self.internalModes.z;
        end
        
        function set.rho(~,~)
            error('This property is readonly.')
        end
        function value = get.rho(self)
            value = self.internalModes.rho;
        end
        
        function set.N2(~,~)
            error('This property is readonly.')
        end
        function value = get.N2(self)
            value = self.internalModes.N2;
        end
        
        function set.rho_z(~,~)
            error('This property is readonly.')
        end
        function value = get.rho_z(self)
            value = self.internalModes.rho_z;
        end
        
        function set.rho_zz(~,~)
            error('This property is readonly.')
        end
        function value = get.rho_zz(self)
            value = self.internalModes.rho_zz;
        end
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Set the normalization and upper boundary condition
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set.normalization(self,norm)
            self.internalModes.normalization = norm;
        end
        function norm = get.normalization(self)
            norm = self.internalModes.normalization;
        end
        
        function set.upperBoundary(self,value)
            self.internalModes.upperBoundary = value;
        end
        function value = get.upperBoundary(self)
            value = self.internalModes.upperBoundary;
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Primary methods to construct the modes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [F,G,h] = ModesAtWavenumber(self, k )
            % Return the normal modes and eigenvalue at a given wavenumber.
            [F,G,h] = self.internalModes.ModesAtWavenumber( k );
        end
        
        function [F,G,h] = ModesAtFrequency(self, omega )
            % Return the normal modes and eigenvalue at a given frequency.
            [F,G,h] = self.internalModes.ModesAtFrequency( omega );
        end
    end
    
    methods (Access = private)
        
        function self = InitTestCase(self, stratification, theMethod, n)
            self.method = theMethod;
            self.isRunningTestCase = 1;
            self.stratification = stratification;
            
            fprintf('InternalModes intialized with %s stratification, %d grid points using %s.\n', self.stratification, n, self.fullMethodName);
            
            lat = 33;
            N0 = 5.2e-3; % reference buoyancy frequency, radians/seconds
            g = 9.81;
            rho0 = 1025;
            if  strcmp(stratification, 'constant')
                self.rhoFunction = @(z) -(N0*N0*rho0/g)*z + rho0;
                self.N2Function = @(z) N0*N0*ones(size(z));
            elseif strcmp(stratification, 'exponential')
                L_gm = 1.3e3; % thermocline exponential scale, meters
                self.rhoFunction = @(z) rho0*(1 + L_gm*N0*N0/(2*g)*(1 - exp(2*z/L_gm)));
                self.N2Function = @(z) N0*N0*exp(2*z/L_gm);
            else
                error('Invalid choice of stratification: you must use constant or exponential');
            end
            zIn = [-5000 0];
            zOut = linspace(zIn(1),0,n)';
            
            if  strcmp(theMethod, 'densitySpectral')
                self.internalModes = InternalModesDensitySpectral(self.rhoFunction,zIn,zOut,lat);
            elseif  strcmp(theMethod, 'wkbSpectral')
                self.internalModes = InternalModesWKBSpectral(self.rhoFunction,zIn,zOut,lat);
            elseif strcmp(theMethod, 'finiteDifference')
                self.internalModes = InternalModesFiniteDifference(self.rhoFunction,zIn,zOut,lat);
            elseif strcmp(theMethod, 'spectral')
                self.internalModes = InternalModesSpectral(self.rhoFunction,zIn,zOut,lat);
            elseif isempty(theMethod)
                self.internalModes = InternalModesWKBSpectral(self.rhoFunction,zIn,zOut,lat);
            else
                error('Invalid method!')
            end
        end
        
        function methodName = fullMethodName( self )
            if  strcmp(self.method, 'densitySpectral')
                methodName = 'Chebyshev polynomials on density coordinates';
            elseif  strcmp(self.method, 'wkbSpectral')
                methodName = 'Chebyshev polynomials on WKB coordinates';
            elseif strcmp(self.method, 'finiteDifference')
                methodName = 'finite differencing';
            elseif strcmp(self.method, 'spectral')
                methodName = 'Chebyshev polynomials';
            elseif isempty(self.method)
                methodName = 'Chebyshev polynomials on WKB coordinates';
            else
                error('Invalid method!')
            end
        end
                
        function self = ShowLowestModesFigure(self,F,G,h)
            figure
            subplot(1,3,1)
            plot(F(:,1:4),self.z, 'LineWidth', 2)
            ylabel('depth (meters)');
            xlabel('(u,v)-modes');
            
            b = subplot(1,3,2);
            plot(G(:,1:4),self.z, 'LineWidth', 2)
            title(b, sprintf('Internal Modes for %s stratification computed using %s\n h = (%.2g, %.2g, %.2g, %.2g)',self.stratification,self.fullMethodName, h(1) , h(2), h(3), h(4) ));
            xlabel('w-modes');
            ytick([]);
            
            subplot(1,3,3)
            plot(sqrt(self.N2),self.z, 'LineWidth', 2), hold on
            if ~isempty(self.N2Function)
                plot(sqrt(self.N2Function(self.z)),self.z, 'LineWidth', 2)
            end
            xlim([0.0 1.1*max(sqrt(self.N2))])
            xlabel('buoyancy frequency');
            ytick([]);
        end
        
        function self = ShowErrorFigure(self, h_error, F_error, G_error, theTitle)
            maxModes = length(h_error);
            
            if strcmp(self.upperBoundary,'free_surface')
                modes = 0:(maxModes-1);
            else
                modes = 1:maxModes;
            end
            
            figure
            plot(modes,h_error, 'g'), ylog
            hold on
            plot(modes,F_error, 'b')
            plot(modes,G_error, 'k')
            xlabel('Mode')
            ylabel('Relative error')
            legend('h', 'F', 'G')
            title(theTitle)
            xlim([0 maxModes])
            ylim([1e-15 1e1])
        end
        
    end
    
end


