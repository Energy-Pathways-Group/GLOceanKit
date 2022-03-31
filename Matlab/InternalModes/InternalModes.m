
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
    %       [F,G,h,omega] = modes.ModesAtWavenumber(0.01);
    %   or frequency,
    %       [F,G,h,k] = modes.ModesAtWavenumber(5*modes.f0);
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
        shouldShowDiagnostics % flag to show diagnostic information, default = 0
        
        latitude % Latitude for which the modes are being computed.
        f0 % Coriolis parameter at the above latitude.
        nModes % Number of modes to be returned.
        
        Lz % Depth of the ocean.
        z % Depth coordinate grid used for all output (same as zOut).
        rho % Density on the z grid.
        N2 % Buoyancy frequency on the z grid, $N^2 = -\frac{g}{\rho(0)} \frac{\partial \rho}{\partial z}$.
        rho_z % First derivative of density on the z grid.
        rho_zz % Second derivative of density on the z grid.
        rho0 % density at the surface (or user specified through constructor args)
        
        lowerBoundary % Lower boundary condition. Either LowerBoundary.freeSlip (default) or LowerBoundary.none.
        upperBoundary % Surface boundary condition. Either UpperBoundary.rigidLid (default) or UpperBoundary.freeSurface.
        normalization % Normalization used for the modes. Either Normalization.(kConstant, omegaConstant, uMax, or wMax).
    end
    
    properties %(Access = private)
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
                
                if isa(rho,'numeric') == true
                   rho = reshape(rho,[],1);
                   zIn = reshape(zIn,[],1);
                   zOut = reshape(zOut,[],1);
                end
                
                isStratificationConstant = InternalModesConstantStratification.IsStratificationConstant(rho,zIn);
                [isStratificationExponential, rho_params] = InternalModesExponentialStratification.IsStratificationExponential(rho,zIn);
                if userSpecifiedMethod == 0 && isStratificationConstant == 1
                    [N0, rho0] = InternalModesConstantStratification.BuoyancyFrequencyFromConstantStratification(rho,zIn);
                    fprintf('Initialization detected that you are using constant stratification with N0=%.1g. The modes will now be computed using the analytical form. If you would like to override this behavior, specify the method parameter.\n', N0);
                    self.internalModes = InternalModesConstantStratification([N0, rho0],zIn,zOut,latitude,extraargs{:});
                elseif userSpecifiedMethod == 0 && isStratificationExponential == 1
                    fprintf('Initialization detected that you are using exponential stratification with N0=%.1g, b=%d. The modes will now be computed using the analytical form. If you would like to override this behavior, specify the method parameter.\n',rho_params(1),round(rho_params(2)));
                    self.internalModes = InternalModesExponentialStratification(rho_params,zIn,zOut,latitude,extraargs{:});
                elseif userSpecifiedMethod == 0
                    % If the user didn't specify a method, try
                    % wkbAdaptiveSpectral first, but if that fails to
                    % produce a reasonable coordinate system, try
                    % z-coordinate spectral.
%                     try
                        self.method = 'wkbAdaptiveSpectral';
                        self.internalModes = InternalModesAdaptiveSpectral(rho,zIn,zOut,latitude,extraargs{:});
%                     catch ME
%                         errorStruct = [];
%                         if ~isempty(ME.cause) && strcmp(ME.cause{1}.identifier,'StretchedGridFromCoordinate:NonMonotonicFunction')
%                             errorStruct = ME.cause{1};
%                         elseif strcmp(ME.identifier,'StretchedGridFromCoordinate:NonMonotonicFunction')
%                             errorStruct = ME;
%                         end
%                         
%                         if ~isempty(errorStruct)
%                             fprintf('%s : %s\n', errorStruct.identifier, errorStruct.message);
%                             fprintf('Switching to InternalModesSpectral.\n');
%                             self.method = 'spectral';
%                             self.internalModes = InternalModesSpectral(rho,zIn,zOut,latitude,extraargs{:});
%                         else
%                             rethrow(ME);
%                         end
%                     end
                elseif  strcmp(self.method, 'densitySpectral')
                    self.internalModes = InternalModesDensitySpectral(rho,zIn,zOut,latitude,extraargs{:});
                elseif  strcmp(self.method, 'wkbSpectral')
                    self.internalModes = InternalModesWKBSpectral(rho,zIn,zOut,latitude,extraargs{:});
                elseif  strcmp(self.method, 'wkbAdaptiveSpectral')
                    self.internalModes = InternalModesAdaptiveSpectral(rho,zIn,zOut,latitude,extraargs{:});
                elseif strcmp(self.method, 'finiteDifference')
                    self.internalModes = InternalModesFiniteDifference(rho,zIn,zOut,latitude,extraargs{:});
                elseif strcmp(self.method, 'spectral')
                    self.internalModes = InternalModesSpectral(rho,zIn,zOut,latitude,extraargs{:});
                elseif strcmp(self.method, 'wkb')
                    self.internalModes = InternalModesWKB(rho,zIn,zOut,latitude,extraargs{:});
                elseif strcmp(self.method, 'wkb-hydrostatic')
                    self.internalModes = InternalModesWKBHydrostatic(rho,zIn,zOut,latitude,extraargs{:});
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
                    theMethod = 'wkbAdaptiveSpectral';
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
            
            % all test cases should just use all modes possible.
            self.internalModes.nModes = length(self.internalModes.z);
            
            [F,G,h] = self.internalModes.ModesAtWavenumber( k );
            
            % y is the true solution, x is the approximated
            errorFunction = @(x,y) max(abs(x-y),[],1)./max(abs(y),[],1);
            
            if strcmp(self.stratification, 'constant')
                imConstant = InternalModesConstantStratification(5.2e-3,[-5000 0],self.z,self.latitude,'nModes',self.nModes);
                imConstant.upperBoundary = self.upperBoundary;
                imConstant.normalization = self.normalization;
                [F_analytical,G_analytical,h_analytical] = imConstant.ModesAtWavenumber( k );
            elseif  strcmp(self.stratification, 'exponential')
                imExponential = InternalModesExponentialStratification([5.2e-3 1300],[-5000 0],self.z,self.latitude,'nModes',self.nModes);
                imExponential.upperBoundary = self.upperBoundary;
                imExponential.normalization = self.normalization;
                [F_analytical,G_analytical,h_analytical] = imExponential.ModesAtWavenumber( k );
            else
                [rhoFunc, ~, zIn] = InternalModes.StratificationProfileWithName(self.stratification);
                imAnalytical = InternalModesAdaptiveSpectral(rhoFunc,zIn,self.z,self.latitude,'nEVP',512,'nModes',self.nModes);
                imAnalytical.upperBoundary = self.upperBoundary;
                imAnalytical.normalization = self.normalization;
                [F_analytical,G_analytical,h_analytical] = imAnalytical.ModesAtWavenumber( k );
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
            
            % all test cases should just use all modes possible.
            self.internalModes.nModes = length(self.internalModes.z);
            
            
            
            % y is the true solution, x is the approximated
            errorFunction = @(x,y) max(max(abs(x-y),[],1)./max(abs(y),[],1),1e-15);
            
            if  strcmp(self.stratification, 'constant')
                imConstant = InternalModesConstantStratification(5.2e-3,[-5000 0],self.z,self.latitude,'nModes',self.nModes);
                imConstant.upperBoundary = self.upperBoundary;
                imConstant.normalization = self.normalization;
                [F_analytical,G_analytical,h_analytical] = imConstant.ModesAtFrequency( omega );
            elseif  strcmp(self.stratification, 'exponential')
                imExponential = InternalModesExponentialStratification([5.2e-3 1300],[-5000 0],self.z,self.latitude,'nModes',self.nModes);
                imExponential.upperBoundary = self.upperBoundary;
                imExponential.normalization = self.normalization;
                [F_analytical,G_analytical,h_analytical] = imExponential.ModesAtFrequency( omega );
            else
                [rhoFunc, ~, zIn] = InternalModes.StratificationProfileWithName(self.stratification);
                imAnalytical = InternalModesAdaptiveSpectral(rhoFunc,zIn,self.z,self.latitude,'nEVP',512,'nModes',self.nModes);
                imAnalytical.upperBoundary = self.upperBoundary;
                imAnalytical.normalization = self.normalization;
                [F_analytical,G_analytical,h_analytical] = imAnalytical.ModesAtFrequency( omega );
            end
            
            self.internalModes.nModes = length(h_analytical);
            [F,G,h] = self.internalModes.ModesAtFrequency( omega );
            
            h_error = errorFunction(h,h_analytical);
            F_error = errorFunction(F,F_analytical);
            G_error = errorFunction(G,G_analytical);
            
            title = sprintf('Relative error of internal modes at omega=%.2gf0 with %d grid points using %s',omega/self.f0,length(self.z),self.fullMethodName);
            self.ShowErrorFigure(h_error,F_error,G_error,title);
        end
              
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Diagnostics
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function set.shouldShowDiagnostics(self,value)
            self.internalModes.shouldShowDiagnostics = value;
        end
        function value = get.shouldShowDiagnostics(self)
            value = self.internalModes.shouldShowDiagnostics;
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
        
        function set.rho0(~,~)
            error('This property is readonly.')
        end
        function value = get.rho0(self)
            value = self.internalModes.rho0;
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
        
        function set.lowerBoundary(self,value)
            self.internalModes.lowerBoundary = value;
        end
        function value = get.lowerBoundary(self)
            value = self.internalModes.lowerBoundary;
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Primary methods to construct the modes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [F,G,h,omega,varargout] = ModesAtWavenumber(self, k, varargin )
            % Return the normal modes and eigenvalue at a given wavenumber.
            if isempty(varargin)
                [F,G,h,omega] = self.internalModes.ModesAtWavenumber( k );
            else
                varargout = cell(size(varargin));
                [F,G,h,omega,varargout{:}] = self.internalModes.ModesAtWavenumber( k, varargin{:} );
            end
        end
        
        function [F,G,h,k,varargout] = ModesAtFrequency(self, omega, varargin )
            % Return the normal modes and eigenvalue at a given frequency.
            if isempty(varargin)
                [F,G,h,k] = self.internalModes.ModesAtFrequency( omega );
            else
                varargout = cell(size(varargin));
                [F,G,h,k,varargout{:}] = self.internalModes.ModesAtFrequency( omega, varargin{:} );
            end
        end
        
        function psi = SurfaceModesAtWavenumber(self, k)
            psi = self.internalModes.SurfaceModesAtWavenumber(k);
        end
        
        function psi = BottomModesAtWavenumber(self, k)
            psi = self.internalModes.BottomModesAtWavenumber(k);
        end
        
        function [m,G] = ProjectOntoGModesAtWavenumber(self, zeta, k)
            
        end
    end
    
    methods (Static)
        
        function N = NumberOfWellConditionedModes(G,varargin)
            % This function can become a rate limiting step if used for
            % each set of returned modes. So a good algorithm is necessary.
            % Otherwise we'd just use,
            %   kappa = InternalModes.ConditionNumberAsFunctionOfModeNumber(G);
            %   N = find(diff(diff(log10(kappa))) > 1e-2,1,'first')+2;
            maxModes = min(size(G));
            minModes = 1;
            
            kappaFull = zeros(maxModes,1);
            kappaFullIsSet = zeros(maxModes,1);
            
            % include the end points, min 4 indices to be helpful
            % maxModes=5, stride = 1
            % maxModes=6, stride = 2
            % maxModes=7, stride = 2
            % maxModes=9, stride = 2
            % maxModes=10, stride = 4
            % maxModes=18, stride = 8
            % maxModes=stride*2+2
            % stride=(maxModes-2)/2
            strideExp = max(0,floor(log2((maxModes-2)/2)));
            lowerBound = minModes;
            upperBound = maxModes;
            while (strideExp >= 0)
                stride = 2^strideExp;
                modeIndices = lowerBound:stride:upperBound;
                if (modeIndices(end) ~= upperBound)
                    modeIndices(end+1) = upperBound;
                end
                
                indicesNeeded = modeIndices(kappaFullIsSet(modeIndices) == 0);
                
                kappaNeeded = InternalModes.ConditionNumberAsFunctionOfModeNumberForModeIndices(G,indicesNeeded);
                kappaFull(indicesNeeded) = kappaNeeded;
                kappaFullIsSet(indicesNeeded) = 1;
                
                kappa = kappaFull(modeIndices);
                
                n = find(diff(diff(log10(kappa))./diff(modeIndices')) > 3e-2,1,'first');
                if isempty(n)
                    %[~,n] = max(diff(diff(log10(kappa))./diff(modeIndices')));
                    n = length(modeIndices)-2;
                    N = modeIndices(n+2);
                else
                    N = modeIndices(n+1);
                end
                lowerBound = modeIndices(n);
                upperBound = modeIndices(n+2);
                strideExp = strideExp-1;
            end
            
        end
        
        function [kappa,modeIndices] = ConditionNumberAsFunctionOfModeNumberForModeIndices(G,modeIndices)
            kappa = zeros(length(modeIndices),1);
            for iIndex=1:length(modeIndices)
               kappa(iIndex) =  cond(G(:,1:modeIndices(iIndex)));
            end
        end
        
        function kappa = ConditionNumberAsFunctionOfModeNumber(G)
            kappa = InternalModes.ConditionNumberAsFunctionOfModeNumberForModeIndices(G,1:min(size(G)));
        end
        
        function [G_tilde, gamma] = RenormalizeForGoodConditioning(G)
            gamma = zeros(1,size(G,2));
            for iMode = 1:size(G,2)
                gamma(iMode) = norm(G(:,iMode));
            end
            gamma = median(gamma)./gamma;
            G_tilde = G .* gamma;
        end
                
        function [rhoFunc, N2Func, zIn] = StratificationProfileWithName(stratification)
            % Returns function handles for several built-in analytical
            % profiles. Options include 'constant', 'exponential',
            % 'pycnocline-constant', 'pycnocline-exponential', and
            % 'latmix-site1'.
            N0 = 5.2e-3; % reference buoyancy frequency, radians/seconds
            g = 9.81;
            rho_0 = 1025;
            if  strcmp(stratification, 'constant')
                rhoFunc = @(z) -(N0*N0*rho_0/g)*z + rho_0;
                N2Func = @(z) N0*N0*ones(size(z));
                zIn = [-5000 0];
            elseif strcmp(stratification, 'exponential')
                L_gm = 1.3e3; % thermocline exponential scale, meters
                rhoFunc = @(z) rho_0*(1 + L_gm*N0*N0/(2*g)*(1 - exp(2*z/L_gm)));
                N2Func = @(z) N0*N0*exp(2*z/L_gm);
                zIn = [-5000 0];
            elseif strcmp(stratification, 'pycnocline-constant')
                % Taken to match the simple test case in Cushman-Roison
                rho_L=1024.63; % these are reversed from the usual definition
                rho_0=1026.45;
                N2_max=2E-3;
                z_ml=-50;
                h = 10;
                L_z = 200;
                
                kappa_p = tanh((L_z-z_ml)/h);
                kappa_m = tanh(-z_ml/h);
                kappa = kappa_p - kappa_m;
                
                beta = (rho_0*(1-h*N2_max*kappa/g) - rho_L)/(L_z - h*kappa);
                gamma = h*(beta-rho_0*N2_max/g);
                delta = rho_L - gamma*kappa_p;
                
                rhoFunc = @(z) delta + beta*(L_z-z)+gamma*tanh((z-z_ml)/h) - L_z;
                N2Func = @(z) -(g/rhoFunc(0))*(-beta + (gamma/h)*sech((z-z_ml)/h).^2);
                
                zIn = [-L_z 0];
            elseif strcmp(stratification, 'pycnocline-exponential')
                delta_p = 75;
                z_s = 200;
                L_s = 500;
                L_d = 2500; %L_d = 1300;
                z_p = -500; %z_p = -1;
                D = 5000;
                
                N0 = 0.5e-3; %N0 = 5.2e-3;
                Np = 10e-3;
                Nd = 3e-4; %Nd = 1.1e-4;
                
                [rhoFunc, N2Func, zIn] = InternalModes.ProfileWithDoubleExponentialPlusPycnocline(rho_0, D, delta_p, z_p, z_s, L_s, L_d, N0, Np, Nd);
            elseif strcmp(stratification, 'latmix-site1')
                delta_p = 275;
                z_s = 200;
                L_s = 500;
                L_d = 2500; %L_d = 1300;
                z_p = -700; %z_p = -1;
                D = 5000;
                
                N0 = 0.5e-3; %N0 = 5.2e-3;
                Np = 5e-3;
                Nd = 2.5e-4; %Nd = 1.1e-4;
                
                rho_0 = 1024.6;
                
                [rho_d, N2_d, zIn] = InternalModes.ProfileWithDoubleExponentialPlusPycnocline(rho_0, D, delta_p, z_p, z_s, L_s, L_d, N0, Np, Nd);
                [rho_s, N2_s, ~] = InternalModes.StratificationProfileWithName('latmix-site1-surface');
                
                rhoFunc = @(z) rho_s(z)-rho_s(0) + rho_d(z);
                N2Func = @(z) N2_s(z) + N2_d(z);
            elseif strcmp(stratification, 'latmix-site1-surface')
                % Recreates the surface mixed layer and associated
                % pycnocline, and then add the deeper thermocline
                delta_p = 0.9;
                L_s = 8;
                L_d = 20;
                z_p = -17;
                z_T = -30;
                D = 5000;
                b = 110;
                
                N0 = 2.8e-3; % Surface
                Nq = 1.4e-2; % Pycnocline portion to exponentials
                Np = 4.7e-2;
                
                [rhoFunc, N2Func, zIn] = InternalModes.ProfileWithDoubleGaussianExponentialPlusPycnocline(rho_0, D, delta_p, z_p, L_s, L_d, z_T, b, N0, Nq, Np);
            elseif strcmp(stratification, 'latmix-site1-exponential')
                % Recreates the surface mixed layer and associated
                % pycnocline, and then decays exponentially below 190m.
                delta_p = 0.9;
                L_s = 8;
                L_d = 100;
                z_p = -17;
                z_T = -190;
                D = 5000;
                b = 1300;
                
                N0 = 2.8e-3; % Surface
                Nq = 1.4e-2; % Pycnocline portion to exponentials
                Np = 4.7e-2;
                
                [rhoFunc, N2Func, zIn] = InternalModes.ProfileWithDoubleGaussianExponentialPlusPycnocline(rho_0, D, delta_p, z_p, L_s, L_d, z_T, b, N0, Nq, Np);
            elseif strcmp(stratification, 'latmix-site1-constant')
                % Recreates the surface mixed layer and associated
                % pycnocline, but then goes constant below 300m.
                delta_p = 0.9;
                L_s = 8;
                L_d = 100;
                z_p = -17; %z_p = -1;
                D = 5000;
                
                N0 = 2.8e-3; % Surface
                Nq = 1.4e-2; % Pycnocline portion to exponentials
                Np = 4.7e-2;
                Nd = 1.75e-3;
                
                [rhoFunc, N2Func, zIn] = InternalModes.ProfileWithDoubleGaussianPlusPycnocline(rho_0, D, delta_p, z_p, L_s, L_d, N0, Nq, Np, Nd);
            else
                error('Invalid choice of stratification. Valid options are: constant, exponential, pycnocline-constant, pycnocline-exponential, or latmix-site1.');
            end
        end
        
        function [rhoFunc, N2Func, zIn] = ProfileWithDoubleExponentialPlusPycnocline(rho_0, D, delta_p, z_p, z_s, L_s, L_d, N0, Np, Nd)            
            A = Np*Np - Nd*Nd*exp(2*(D+z_p)/L_d);
            B = (Nd*Nd*exp(2*(D+z_p)/L_d) - N0*N0)/(exp(-2*(z_p-z_s)/L_s) - exp(2*z_s/L_s));
            C = N0*N0-B*exp(2*z_s/L_s);
            E = Nd*Nd*exp(2*(D+z_p)/L_d);
            
            g = 9.81;
            rho_surface = @(z) rho_0*(1 - (L_s*B/(2*g)) * (exp(2*z_s/L_s) - exp(-2*(z-z_s)/L_s)) - C*z/g);
            rho_deep = @(z) rho_surface(z_p) + (rho_0*L_d*E/(2*g))*(1-exp(2*(z-z_p)/L_d));
            rho_p = @(z) -(A*rho_0*delta_p/g)*(tanh( (z-z_p)/delta_p) - 1);
            
            rhoFunc = @(z) (z>=z_p).*rho_surface(z) + (z<z_p).*rho_deep(z) + rho_p(z);
            
            N2_surface = @(z) B*exp(-2*(z-z_s)/L_s) + C;
            N2_deep = @(z) E*exp(2*(z-z_p)/L_d);
            N2_p = @(z) A*sech( (z-z_p)/delta_p ).^2;
            
            N2Func = @(z) (z>=z_p).*N2_surface(z) + (z<z_p).*N2_deep(z) + N2_p(z);
            zIn = [-D 0];
        end
        
        function [rhoFunc, N2Func, zIn] = ProfileWithDoubleGaussianPlusPycnocline(rho_0, D, delta_p, z_p, L_s, L_d, N0, Nq, Np, Nd)
            A = Np*Np - Nq*Nq;
            B = (Nq*Nq - N0*N0)/(1 - exp(-2*z_p^2/L_s^2));
            C = N0*N0-B*exp(-2*z_p^2/L_s);
            E = (Nq*Nq - Nd*Nd)/( 1 - exp(-2*(-D-z_p)^2/L_d^2) );
            F = Nd*Nd - E*exp( -2*(-D-z_p)^2/L_d^2 );
            
            g = 9.81;
            rho_surface = @(z) rho_0*(1 - (L_s*B/(2*g)) * sqrt(pi/2) *( erf(sqrt(2)*z_p/L_s) - erf(-sqrt(2)*(z-z_p)/L_s) ) - C*z/g);
            rho_deep = @(z) rho_surface(z_p) - (rho_0*L_d*E/(2*g)) * sqrt(pi/2) * (erf(sqrt(2)*(z-z_p)/L_d)) - rho_0*F*(z-z_p)/g;
            rho_p = @(z) -(A*rho_0*delta_p/g)*(tanh( (z-z_p)/delta_p) - 1);
            
            rhoFunc = @(z) (z>=z_p).*rho_surface(max(z,z_p)) + (z<z_p).*rho_deep(z) + rho_p(z);
            
            N2_surface = @(z) B*exp(-2*(z-z_p).^2/L_s^2) + C;
            N2_deep = @(z) E*exp(-2*(z-z_p).^2/L_d^2) + F;
            N2_p = @(z) A*sech( (z-z_p)/delta_p ).^2;
            
            N2Func = @(z) (z>=z_p).*N2_surface(z) + (z<z_p).*N2_deep(z) + N2_p(z);
            zIn = [-D 0];
        end
        
        function [rhoFunc, N2Func, zIn] = ProfileWithDoubleGaussianExponentialPlusPycnocline(rho_0, D, delta_p, z_p, L_s, L_d, z_T, L_deep, N0, Nq, Np)
            A = Np*Np - Nq*Nq;
            B = (Nq*Nq - N0*N0)/(1 - exp(-2*z_p^2/L_s^2));
            C = N0*N0-B*exp(-2*z_p^2/L_s);
            alpha = -(2*L_deep*(z_T-z_p)/(L_d^2))*exp(-2*(z_T-z_p)^2/(L_d^2) - 2*z_T/L_deep);
            gamma = alpha*exp(2*z_T/L_deep)-exp(-2*(z_T-z_p)^2/(L_d^2));
            E = Nq*Nq/( 1 + gamma );
            F = E*gamma;
            G = E*alpha;
            
            g = 9.81;
            rho_surface = @(z) rho_0*(1 - (L_s*B/(2*g)) * sqrt(pi/2) *( erf(sqrt(2)*z_p/L_s) - erf(-sqrt(2)*(z-z_p)/L_s) ) - C*z/g);
            rho_mid = @(z) rho_surface(z_p) - (rho_0*L_d*E/(2*g)) * sqrt(pi/2) * (erf(sqrt(2)*(z-z_p)/L_d)) - rho_0*F*(z-z_p)/g;
            rho_deep = @(z) rho_mid(z_T) - (rho_0*L_deep/(2*g))*G*(exp(2*z/L_deep)-exp(2*z_T/L_deep));
            rho_p = @(z) -(A*rho_0*delta_p/g)*(tanh( (z-z_p)/delta_p) - 1);
            
            rhoFunc = @(z) (z>=z_p).*rho_surface(max(z,z_p)) + (z<z_p & z > z_T).*rho_mid(z) + (z<=z_T).*rho_deep(z) + rho_p(z);
            
            N2_surface = @(z) B*exp(-2*(z-z_p).^2/L_s^2) + C;
            N2_mid = @(z) E*exp(-2*(z-z_p).^2/L_d^2) + F;
            N2_deep = @(z) G*exp(2*z/L_deep);
            N2_p = @(z) A*sech( (z-z_p)/delta_p ).^2;
            
            N2Func = @(z) (z>=z_p).*N2_surface(z) + (z<z_p & z >z_T).*N2_mid(z) + (z<=z_T).*N2_deep(z) + N2_p(z);
            zIn = [-D 0];
        end
    end
    
    methods (Access = private)
        
        function self = InitTestCase(self, stratification, theMethod, n)
            self.method = theMethod;
            self.isRunningTestCase = 1;
            self.stratification = stratification;
            
            fprintf('InternalModes intialized with %s stratification, %d grid points using %s.\n', self.stratification, n, self.fullMethodName);
            
            lat = 33;
            
            [self.rhoFunction, self.N2Function, zIn] = InternalModes.StratificationProfileWithName(stratification);
            zOut = linspace(min(zIn),max(zIn),n)';

            if  strcmp(theMethod, 'densitySpectral')
                self.internalModes = InternalModesDensitySpectral(self.rhoFunction,zIn,zOut,lat,'nEVP', n);
            elseif  strcmp(theMethod, 'wkbSpectral')
                self.internalModes = InternalModesWKBSpectral(self.rhoFunction,zIn,zOut,lat,'nEVP', n);
            elseif  strcmp(theMethod, 'wkbAdaptiveSpectral')
                self.internalModes = InternalModesAdaptiveSpectral(self.rhoFunction,zIn,zOut,lat,'nEVP', n);
            elseif strcmp(theMethod, 'finiteDifference')
                self.internalModes = InternalModesFiniteDifference(self.rhoFunction,zIn,zOut,lat);
            elseif strcmp(theMethod, 'spectral')
                self.internalModes = InternalModesSpectral(self.rhoFunction,zIn,zOut,lat,'nEVP', n);
            elseif strcmp(self.method, 'wkb')
                self.internalModes = InternalModesWKB(self.rhoFunction,zIn,zOut,lat);
            elseif strcmp(self.method, 'wkb-hydrostatic')
                self.internalModes = InternalModesWKBHydrostatic(self.rhoFunction,zIn,zOut,lat);
            elseif isempty(theMethod)
                self.internalModes = InternalModesWKBSpectral(self.rhoFunction,zIn,zOut,lat,'nEVP', n);
            else
                error('Invalid method!')
            end
        end
        
        function methodName = fullMethodName( self )
            if  strcmp(self.method, 'densitySpectral')
                methodName = 'Chebyshev polynomials on density coordinates';
            elseif  strcmp(self.method, 'wkbSpectral')
                methodName = 'Chebyshev polynomials on WKB coordinates';
            elseif  strcmp(self.method, 'wkbAdaptiveSpectral')
                methodName = 'Chebyshev polynomials on WKB adaptive coordinates';
            elseif strcmp(self.method, 'finiteDifference')
                methodName = 'finite differencing';
            elseif strcmp(self.method, 'spectral')
                methodName = 'Chebyshev polynomials';
            elseif strcmp(self.method, 'wkb')
                methodName = 'WKB approximation';
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
            yticks([]);
            
            subplot(1,3,3)
%             if any(self.N2 < 0)
                plot(self.N2,self.z, 'LineWidth', 2), hold on
                xlim([0.9*min(self.N2) 1.1*max(self.N2)])
                xlabel('N^2');
%             else
%                 plot(sqrt(self.N2),self.z, 'LineWidth', 2), hold on
%                 if ~isempty(self.N2Function)
%                     plot(sqrt(self.N2Function(self.z)),self.z, 'LineWidth', 2)
%                 end
%                 xlim([0.0 1.1*max(sqrt(self.N2))])
%                 xlabel('buoyancy frequency');
%             end
            
            yticks([]);
        end
        
        function self = ShowErrorFigure(self, h_error, F_error, G_error, theTitle)
            maxModes = length(h_error);
            
            if self.upperBoundary == UpperBoundary.freeSurface
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


