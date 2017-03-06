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
classdef InternalModes < handle
    properties (Access = public)
        method
        internalModes % instance of actual internal modes class that's doing all the work
    end
    
    properties (Dependent) 
        latitude
        f0
        
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
            self.method = 'stretchedSpectral';
            
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
                elseif  strcmp(self.method, 'stretchedSpectral')
                    self.internalModes = InternalModesStretchedSpectral(rho,zIn,zOut,latitude,extraargs{:});
                elseif strcmp(self.method, 'finiteDifference')
                    self.internalModes = InternalModesFiniteDifference(rho,zIn,zOut,latitude,extraargs{:});
                elseif strcmp(self.method, 'spectral')
                    self.internalModes = InternalModesSpectral(rho,zIn,zOut,latitude,extraargs{:});
                else
                    error('Invalid method!')
                end
            else
                fprintf('InternalModes intialized for test case usage.\n');
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function self = ShowTestCase(self, stratification, theMethod )
            lat = 33;
            N0 = 5.2e-3; % reference buoyancy frequency, radians/seconds
            g = 9.81;
            rho0 = 1025;
            if  strcmp(stratification, 'constant')
                rhoFunction = @(z) -(N0*N0*rho0/g)*z + rho0;
                N2Function = @(z) N0*N0*ones(size(z));
            elseif strcmp(stratification, 'exponential')
                L_gm = 1.3e3; % thermocline exponential scale, meters
                rhoFunction = @(z) rho0*(1 + L_gm*N0*N0/(2*g)*(1 - exp(2*z/L_gm)));
                N2Function = @(z) N0*N0*exp(2*z/L_gm);
            else
                error('Invalid choice of stratification: you must use constant or exponential');
            end
            n = 64;
            zIn = [-5000 0];
            zOut = linspace(zIn(1),0,n)';
             
            if  strcmp(theMethod, 'densitySpectral')
                self.internalModes = InternalModesDensitySpectral(rhoFunction,zIn,zOut,lat);
                methodName = 'Chebyshev polynomials on density coordinates';
            elseif  strcmp(theMethod, 'stretchedSpectral')
                self.internalModes = InternalModesStretchedSpectral(rhoFunction,zIn,zOut,lat);
                methodName = 'Chebyshev polynomials on WKB coordinates';
            elseif strcmp(theMethod, 'finiteDifference')
                self.internalModes = InternalModesFiniteDifference(rhoFunction,zIn,zOut,lat);
                methodName = 'finite differencing';
            elseif strcmp(theMethod, 'spectral')
                self.internalModes = InternalModesSpectral(rhoFunction,zIn,zOut,lat);
                methodName = 'Chebyshev polynomials';
            else
                self.internalModes = InternalModesStretchedSpectral(rhoFunction,zIn,zOut,lat);
                methodName = 'Chebyshev polynomials on WKB coordinates';
            end
            
            k=0.0;
            [F,G,h] = self.internalModes.ModesAtWavenumber( k );
            
            figure
            subplot(1,3,1)
            plot(F(:,1:4),self.z, 'LineWidth', 2)
            ylabel('depth (meters)');
            xlabel('(u,v)-modes');
            
            b = subplot(1,3,2);
            plot(G(:,1:4),self.z, 'LineWidth', 2)
            title(b, sprintf('Internal Modes for %s stratification computed using %s\n h = (%.2g, %.2g, %.2g, %.2g)',stratification,methodName, h(1) , h(2), h(3), h(4) ));
            xlabel('w-modes');
            ytick([]);
            
            subplot(1,3,3)
            plot(self.N2,self.z, 'LineWidth', 2), hold on
            plot(N2Function(self.z),self.z, 'LineWidth', 2)
            xlim([0.0 1.1*max(self.N2)])
            xlabel('buoyancy frequency');
            ytick([]);
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Vertical grid and derivatives of the density profile on that grid
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
end

