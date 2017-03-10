%
%
%   modes = InternalModes(rho,z,z_out,latitude)
%   'rho' must be a vector with length matching 'z'. For finite
%   differencing, z_out is determined with interpolation, with spectral
%   methods, it is projected onto the output grid.
%
%   modes = InternalModes(rho,z_domain,z_out,latitude)
%   'rho' must be a function handle. z_domain must be an array with two
%   values: z_domain = [z_bottom z_surface];
%
%   modes = InternalModes('constant','wkbSpectral',128)
%   initializes with given stratification, method, and grid points
%
classdef InternalModes < handle
    properties (Access = public)
        method
        internalModes % instance of actual internal modes class that's doing all the work
    end
    
    properties (Access = private)
       isRunningTestCase = 0;
       stratification = 'user specified';
       rhoFunction
       N2Function
    end
    
    properties (Dependent) 
        latitude
        f0
        nModes
        
        Lz
        z
        rho
        N2
        rho_z
        rho_zz
        
        upperBoundary
        normalization
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function self = InternalModes(varargin)    
            self.method = 'wkbSpectral';
            
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
                
                if  strcmp(self.method, 'densitySpectral')
                    self.internalModes = InternalModesDensitySpectral(rho,zIn,zOut,latitude,extraargs{:});
                elseif  strcmp(self.method, 'wkbSpectral')
                    self.internalModes = InternalModesWKBSpectral(rho,zIn,zOut,latitude,extraargs{:});
                elseif strcmp(self.method, 'finiteDifference')
                    self.internalModes = InternalModesFiniteDifference(rho,zIn,zOut,latitude,extraargs{:});
                elseif strcmp(self.method, 'spectral')
                    self.internalModes = InternalModesSpectral(rho,zIn,zOut,latitude,extraargs{:});
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Useful functions
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function self = ShowLowestModesAtWavenumber( self, k )
            [F,G,h] = self.internalModes.ModesAtWavenumber( k );
            self.ShowLowestModesFigure(F,G,h);
        end
        
        function self = ShowLowestModesAtFrequency( self, omega )
            [F,G,h] = self.internalModes.ModesAtFrequency( omega );
            self.ShowLowestModesFigure(F,G,h);
        end
        
        function self = ShowRelativeErrorAtWavenumber( self, k )
            if self.isRunningTestCase == 0
                error('Cannot show relative error for user specified stratification.\n');
            end
            
            [F,G,h] = self.internalModes.ModesAtWavenumber( k );
            
            % y is the true solution, x is the approximated
            errorFunction = @(x,y) max(abs(x-y),[],1)./max(abs(y),[],1);
            
            if  strcmp(self.stratification, 'constant')
                [F_analytical,G_analytical,h_analytical] = self.ConstantStratificationModesAtWavenumber(k);
            elseif strcmp(self.stratification, 'exponential')
                imAnalytical = InternalModesWKBSpectral(self.rhoFunction,[-5000 0],self.z,self.latitude,'nEVP',512);
                imAnalytical.nModes = maxModes;
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
        
        function self = ShowRelativeErrorAtFrequency( self, omega )
            if self.isRunningTestCase == 0
                error('Cannot show relative error for user specified stratification.\n');
            end
            
            [F,G,h] = self.internalModes.ModesAtFrequency( omega );
            
            % y is the true solution, x is the approximated
            errorFunction = @(x,y) max(abs(x-y),[],1)./max(abs(y),[],1);
            
            if  strcmp(self.stratification, 'constant')
                [F_analytical,G_analytical,h_analytical] = self.ConstantStratificationModesAtFrequency( omega );
            elseif strcmp(self.stratification, 'exponential')
                imAnalytical = InternalModesWKBSpectral(self.rhoFunction,[-5000 0],self.z,self.latitude,'nEVP',512);
                imAnalytical.nModes = maxModes;
                [F_analytical,G_analytical,h_analytical] = imAnalytical.ModesAtFrequency( k );
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
            [F,G,h] = self.internalModes.ModesAtWavenumber( k );
        end
        
        function [F,G,h] = ModesAtFrequency(self, omega )
            [F,G,h] = self.internalModes.ModesAtFrequency( omega );
        end
    end
    
    methods (Access = private)
        
        function [F,G,h] = ConstantStratificationModesAtWavenumber(self, k)
            N0 = 5.2e-3; g = 9.81;
            k_z = (1:self.nModes)*pi/self.Lz;
            h = (N0*N0 - self.f0*self.f0)./(g*(k*k+k_z.*k_z));
            [F,G] = self.ConstantStratificationModesWithEigenvalue(k_z,h);
        end
        
        function [F,G,h] = ConstantStratificationModesAtFrequency(self, omega)
            N0 = 5.2e-3; g = 9.81;
            k_z = (1:self.nModes)*pi/self.Lz;
            h = (N0*N0 - omega*omega)./(g * k_z.*k_z);
            [F,G] = self.ConstantStratificationModesWithEigenvalue(k_z,h);
        end
        
        % k_z and h should be of size [1, nModes]
        % [F,G] will return with size [length(z), nModes]
        function [F,G] = ConstantStratificationModesWithEigenvalue(self, k_z, h)
            N0 = 5.2e-3; % reference buoyancy frequency, radians/seconds
            g = 9.81;
            if strcmp(self.normalization, 'const_G_norm')
                G = sqrt(2*g/(self.Lz*(N0*N0-self.f0*self.f0))) * sin(k_z .* self.z);
                F = sqrt(2*g/(self.Lz*(N0*N0-self.f0*self.f0))) * repmat(h.*k_z,length(self.z),1) .* cos(k_z .* self.z);
            elseif strcmp(self.normalization, 'const_F_norm')
                G = sqrt(2) * sin(k_z.*self.z) ./ repmat(h.*k_z,length(self.z),1);
                F = sqrt(2) * cos(k_z.*self.z);
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
            plot(self.N2,self.z, 'LineWidth', 2), hold on
            if ~isempty(self.N2Function)
                plot(self.N2Function(self.z),self.z, 'LineWidth', 2)
            end
            xlim([0.0 1.1*max(self.N2)])
            xlabel('buoyancy frequency');
            ytick([]);
        end
        
        function self = ShowErrorFigure(self, h_error, F_error, G_error, theTitle)
            maxModes = length(h_error);
            
            figure
            plot(h_error, 'g'), ylog
            hold on
            plot(F_error, 'b')
            plot(G_error, 'k')
            xlabel('Mode')
            ylabel('Relative error')
            legend('h', 'F', 'G')
            title(theTitle)
            xlim([0 maxModes])
            ylim([1e-15 1e1])
        end
        
    end
    
end


