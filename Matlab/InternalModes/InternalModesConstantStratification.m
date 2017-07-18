classdef InternalModesConstantStratification < InternalModesBase
    properties (Access = public)
        N0
        rho
        N2
        rho_z
        rho_zz
    end
    
    methods
        function self = InternalModesConstantStratification(rho, z_in, z_out, latitude, varargin) 
            self@InternalModesBase(rho,z_in,z_out,latitude, varargin{:});
            self.N0 = InternalModesConstantStratification.BuoyancyFrequencyFromConstantStratification(rho,z_in,self.rho0,self.g);
            
            rhoFunction = @(z) -(self.N0*self.N0*self.rho0/self.g)*z + self.rho0;
            N2Function = @(z) self.N0*self.N0*ones(size(z));
            self.rho = rhoFunction(self.z);
            self.N2 = N2Function(self.z);
            self.rho_z = -(self.N0*self.N0*self.rho0/self.g)*ones(size(self.z));
            self.rho_zz = zeros(size(self.z));
            
            fprintf('Using the analytical form for constant stratification N0=%.7g\n',self.N0);
        end
                
        function [F,G,h] = ModesAtWavenumber(self, k )
            k_z = (1:self.nModes)*pi/self.Lz;
            h = (self.N0*self.N0 - self.f0*self.f0)./(self.g*(k*k+k_z.*k_z));
            

            if strcmp(self.upperBoundary,'free_surface')
               k_z(2:end) = k_z(1:end-1); h(2:end) = h(1:end-1);
               k_z(1) = 0; h(1) = self.Lz;
               
               m = zeros(1,self.nModes);
               h2 = zeros(1,self.nModes);
               
               f = @(m) m*(self.N0*self.N0 - self.f0*self.f0)./(self.g*(k*k-m*m)) - tanh(m*self.Lz);
               m(1) = fzero(f, max(k-1/self.Lz,1/self.Lz));
               h2(1) = (self.N0*self.N0 - self.f0*self.f0)./(self.g*(k*k - m(1)*m(1) ));
               fprintf('For hyperbolic: (m,h) = (%.2g, %.2g)\n', m(1), h2(1));
               
               f = @(m) m*(self.N0*self.N0 - self.f0*self.f0)./(self.g*(k*k+m*m)) - tan(m*self.Lz);
               m(1) = fzero(f, max(k-1/self.Lz,1/self.Lz));
               h2(1) = (self.N0*self.N0 - self.f0*self.f0)./(self.g*(k*k + m(1)*m(1) ));
               fprintf('For trig: (m,h) = (%.2g, %.2g)\n', m(1), h2(1));
               
               f = @(m) m*(self.N0*self.N0 - self.f0*self.f0)./(self.g*(k*k+m*m)) - tan(m*self.Lz);
               x0 = (pi/self.Lz)*(1:self.nModes);
               for i=1:self.nModes-1
                   m(i+1) = fzero(f,x0(i));
               end
               h2(2:end) = (self.N0*self.N0 - self.f0*self.f0)./(self.g*(k*k - m(2:end).*m(2:end) ));
            end
            [F,G] = self.ConstantStratificationModesWithEigenvalue(k_z,h);
        end
        
        function [F,G,h] = ModesAtFrequency(self, omega )
            k_z = (1:self.nModes)*pi/self.Lz;
            h = (self.N0*self.N0 - omega*omega)./(self.g * k_z.*k_z);
            if strcmp(self.upperBoundary,'free_surface')
                k_z(2:end) = k_z(1:end-1); h(2:end) = h(1:end-1);
                k_z(1) = 0; h(1) = self.Lz;
            end
            [F,G] = self.ConstantStratificationModesWithEigenvalue(k_z,h);
        end
        
        % k_z and h should be of size [1, nModes]
        % [F,G] will return with size [length(z), nModes]
        function [F,G] = ConstantStratificationModesWithEigenvalue(self, k_z, h)
            N0_ = self.N0; % reference buoyancy frequency, radians/seconds
            g_ = self.g;
            if strcmp(self.normalization, 'const_G_norm')
                G = sqrt(2*g_/(self.Lz*(N0_*N0_-self.f0*self.f0))) * sin(k_z .* self.z);
                F = sqrt(2*g_/(self.Lz*(N0_*N0_-self.f0*self.f0))) * repmat(h.*k_z,length(self.z),1) .* cos(k_z .* self.z);
                if strcmp(self.upperBoundary,'free_surface')
                    A = sqrt(3*g_/( (N0_*N0_ - self.f0*self.f0)*self.Lz*self.Lz*self.Lz));
                    G(:,1) = A*(self.z + self.Lz);
                    F(:,1) = A*self.Lz*ones(size(self.z));
                end
            elseif strcmp(self.normalization, 'const_F_norm')
                G = sqrt(2) * sin(k_z.*self.z) ./ repmat(h.*k_z,length(self.z),1);
                F = sqrt(2) * cos(k_z.*self.z);
                if strcmp(self.upperBoundary,'free_surface')
                    A = sqrt(1/(self.Lz*self.Lz));
                    G(:,1) = A*(self.z + self.Lz);
                    F(:,1) = A*self.Lz*ones(size(self.z));
                end
            end
        end
        
    end
    
    methods (Access = protected)
        function self = InitializeWithGrid(self, rho, zIn)
            if isempty(self.nModes) || self.nModes < 1
                self.nModes = floor(length(self.z));
            end
        end
        
        function self = InitializeWithFunction(self, rho, zMin, zMax, zOut)
            if isempty(self.nModes) || self.nModes < 1
                self.nModes = floor(length(self.z));
            end
        end
    end
    
    methods (Static)
        function flag = IsStratificationConstant(rho,z_in)
            if isa(rho,'function_handle') == true
                if numel(z_in) ~= 2
                    error('When using a function handle, z_domain must be an array with two values: z_domain = [z_bottom z_surface];')
                end
                rho0 = rho(max(z_in));
                max_ddrho = max(abs(diff(diff(rho(linspace(min(z_in),max(z_in),100))/rho0 ))));
            elseif isa(rho,'numeric') == true
                rho0 = max(rho);
                max_ddrho = max(abs(diff(diff( rho/rho0 ))));
            end
            
            flag = max_ddrho < 1e-7;
        end
        
        function N0 = BuoyancyFrequencyFromConstantStratification(rho,z_in,rho0,g)
            if isa(rho,'function_handle') == true
                drhodz = (rho(max(z_in)) - rho(min(z_in)))/( max(z_in) - min(z_in) );
            elseif isa(rho,'numeric') == true
                drhodz = (min(rho)-max(rho))/(max(z_in)-min(z_in));
            end
            N0 = sqrt(-g*drhodz/rho0);
        end
        
    end
end