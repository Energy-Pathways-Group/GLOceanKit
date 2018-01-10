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
                
        function [F,G,h,omega] = ModesAtWavenumber(self, k )
            k_z = (1:self.nModes)*pi/self.Lz;
            if self.upperBoundary == UpperBoundary.freeSurface % add the free surface correction to the vertical wavenumber
                for i=1:self.nModes
                    f = @(xi) (xi+i*pi)*(self.N0*self.N0 - self.f0*self.f0)*self.Lz - self.g*(k*k*self.Lz*self.Lz+(xi+i*pi).*(xi+i*pi)).*tan(xi);
                    k_z(i) = k_z(i) + fzero(f,0)/self.Lz;
                end
            end
            h = (self.N0*self.N0 - self.f0*self.f0)./(self.g*(k*k+k_z.*k_z)); 
            
            % Now compute the baroclinic modes
            [F,G] = self.BaroclinicModesWithEigenvalue(k_z,h);
            
            if self.upperBoundary == UpperBoundary.freeSurface
                % Make some room for the barotropic mode
                h = circshift(h,1);
                F = circshift(F,1,2);
                G = circshift(G,1,2);
                
                [F0,G0,h0] = self.BarotropicModeAtWavenumber(k);
                h(1) = h0;
                F(:,1) = F0;
                G(:,1) = G0;
            end
            omega = self.omegaFromK(h,k);
        end
        
        function [F,G,h,k] = ModesAtFrequency(self, omega )
            k_z = (1:self.nModes)*pi/self.Lz;
            if self.upperBoundary == UpperBoundary.freeSurface % add the free surface correction to the vertical wavenumber
                for i=1:self.nModes
                    f = @(xi) self.g*tan(xi)/(xi+i*pi) - (self.N0*self.N0 - omega*omega)*self.Lz/((xi+i*pi)^2);
                    k_z(i) = k_z(i) + fzero(f,0)/self.Lz;
                end
            end
            h = (self.N0*self.N0 - omega*omega)./(self.g * k_z.*k_z);
            if (omega >= self.N0)
               k_z = zeros(size(k_z));
               h = zeros(size(h));
            end
            
            [F,G] = self.BaroclinicModesWithEigenvalue(k_z,h);
            
            if self.upperBoundary == UpperBoundary.freeSurface
                % Make some room for the barotropic mode
                h = circshift(h,1);
                F = circshift(F,1,2);
                G = circshift(G,1,2);
                
                [F0,G0,h0] = self.BarotropicModeAtFrequency(omega);
                h(1) = h0;
                F(:,1) = F0;
                G(:,1) = G0;
            end
            k = self.kFromOmega(h,omega);
        end
        
        function [psi] = SurfaceModesAtWavenumber(self, k)
            % size(psi) = [size(k); length(z)]
            sizeK = size(k);
            if length(sizeK) == 2 && sizeK(2) == 1
                sizeK(2) = [];
            end
            k = reshape(k,[],1);
            zCoord = reshape(self.z,1,[]);
            
            lambda = k * (self.N0/self.f0);
            % dividing by sinh is not well conditioned.
            % psi = (1/self.N0) * cosh( (zCoord + self.Lz) .* lambda ) ./ (k .* sinh(lambda*self.Lz));
            psi = (1/self.N0) * (exp( zCoord .* lambda ) + exp( -(zCoord + 2*self.Lz) .* lambda ) ) ./ (k .* (1-exp(-2*lambda*self.Lz)));
            
            sizeK(end+1) = length(self.z);
            psi = reshape(psi,sizeK);
        end
        
        function [psi] = BottomModesAtWavenumber(self, k)
            % size(psi) = [size(k); length(z)]
            sizeK = size(k);
            if length(sizeK) == 2 && sizeK(2) == 1
                sizeK(2) = [];
            end
            k = reshape(k,[],1);
            zCoord = reshape(self.z,1,[]);
            
            lambda = k * (self.N0/self.f0);
            % This form (from the manuscript) is not well conditioned,
            % so we rewrite it in a form that is.
            % psi = (1/self.N0) * ( sinh( (zCoord + self.Lz) .* lambda ) -  cosh( (zCoord + self.Lz) .* lambda ) ./ tanh(lambda*self.Lz)) ./k;
            A = -exp( -lambda*self.Lz )./(1-exp( -2*lambda*self.Lz ));
            B = -exp( -2*lambda*self.Lz )./(1-exp( -2*lambda*self.Lz )) - 1;
            psi = (A.*exp( lambda*zCoord ) + B.*exp( -lambda*(zCoord + self.Lz) ))./(self.N0 * k);
           
            sizeK(end+1) = length(self.z);
            psi = reshape(psi,sizeK);
        end
        
        % k_z and h should be of size [1, nModes]
        % [F,G] will return with size [length(z), nModes]
        function [F,G] = BaroclinicModesWithEigenvalue(self, k_z, h)
            N0_ = self.N0; % reference buoyancy frequency, radians/seconds
            g_ = self.g;
            j = 1:self.nModes;
            switch self.normalization
                case Normalization.kConstant
                    A = (-1).^j .* sqrt(g_./((N0_*N0_-self.f0*self.f0) .* (self.Lz/2 - sin(2*k_z*self.Lz)./(4*k_z))));
                    A = (-1).^j .* (sin(k_z*self.Lz).^2 + (self.Lz/(2*self.g))*(self.N0*self.N0 - self.f0*self.f0)*(1-sin(2*k_z*self.Lz)./(2*k_z*self.Lz))).^(-1/2);
                case Normalization.omegaConstant
                    A = (-1).^j./( h .* k_z .* sqrt(1/2 + sin(2*k_z*self.Lz)./(4*k_z*self.Lz)));
                case Normalization.wMax
                    A = (-1).^j;
                case Normalization.uMax
                    A = (-1).^j./(h.*k_z);
            end
            G = A .*  sin(k_z .* (self.z + self.Lz));
            F = A .*  repmat(h.*k_z,length(self.z),1) .* cos(k_z .* (self.z + self.Lz));
        end
        
        function [F0,G0,h0] = BarotropicModeAtWavenumber(self, k)
            k_star = sqrt( (self.N0*self.N0 - self.f0*self.f0)/(self.g*self.Lz) );
                
            if (abs(k-k_star)/k_star < 1e-6) % transition (linear) solution
                solutionType = 'linear';
                h0 = self.Lz;
                k_z = 0;
            elseif k > k_star % hyperbolic solution
                solutionType = 'hyperbolic';
                f = @(q) self.Lz*(self.N0*self.N0 - self.f0*self.f0) - (1./q).*(self.g*(k*k*self.Lz*self.Lz-q.*q)).*tanh(q);
                k_initial = sqrt( k*k*self.Lz*self.Lz - self.Lz*(self.N0*self.N0 - self.f0*self.f0)/self.g);
                k_z = fzero(f, k_initial)/self.Lz;
                h0 = (self.N0*self.N0 - self.f0*self.f0)./(self.g*(k*k - k_z*k_z ));
            elseif k < k_star % trig solution
                solutionType = 'trig';
                f = @(q) self.Lz*(self.N0*self.N0 - self.f0*self.f0) - (1./q).*(self.g*(k*k*self.Lz*self.Lz+q.*q)).*tan(q);
                k_initial = sqrt( - k*k*self.Lz*self.Lz + self.Lz*(self.N0*self.N0 - self.f0*self.f0)/self.g);
                k_z = fzero(f, k_initial)/self.Lz;
                h0 = (self.N0*self.N0 - self.f0*self.f0)./(self.g*(k*k + k_z*k_z ));        
            end
            
            [F0,G0] = self.BarotropicMode(solutionType, k_z, h0);
        end
        
        function [F0,G0,h0] = BarotropicModeAtFrequency(self, omega)
            if (abs(omega-self.N0)/self.N0 < 1e-6)
                solutionType = 'linear';
                h0 = self.Lz;
                k_z = 0;
            elseif omega > self.N0 % hyperbolic solution
                solutionType = 'hyperbolic';
                f = @(q) self.Lz*(omega*omega - self.N0*self.N0) - self.g * q .* tanh(q);
                k_z = fzero(f, sqrt(self.Lz*(omega*omega - self.N0*self.N0)/self.g) )/self.Lz;
                h0 = (omega*omega - self.N0*self.N0)./(self.g*k_z*k_z);
            elseif omega < self.N0 % trig solution
                solutionType = 'trig';
                f = @(q) self.Lz*(self.N0*self.N0 - omega*omega) - self.g * q .* tan(q);
                k_z = fzero(f, sqrt(self.Lz*(self.N0*self.N0 - omega*omega)/self.g) )/self.Lz;
                h0 = (self.N0*self.N0 - omega*omega)./(self.g * k_z.*k_z);
            end
            
            [F0,G0] = self.BarotropicMode(solutionType, k_z, h0);
        end
        
        function [F0,G0] = BarotropicMode(self, solutionType, k_z, h0)
            % It's safer to do a switch on solutionType, rather than check
            % that omega *or* k are equal to N0, k_star within tolerance.
            if strcmp(solutionType, 'linear')
                switch self.normalization
                    case Normalization.kConstant
                        A = sqrt(3*self.g/( (self.N0*self.N0 - self.f0*self.f0)*self.Lz*self.Lz*self.Lz));
                        A = 1/(self.Lz * sqrt(1 + (self.N0*self.N0 - self.f0*self.f0)*self.Lz/(2*self.g)));
                    case Normalization.omegaConstant
                        A = 1/self.Lz;
                    case Normalization.wMax
                        A = 1/self.Lz;
                    case Normalization.uMax
                        A = 1/self.Lz;
                end
                G0 = A*(self.z + self.Lz);
                F0 = A*self.Lz*ones(size(self.z));
            elseif strcmp(solutionType, 'hyperbolic')
                switch self.normalization
                    case Normalization.kConstant
                        A = sqrt( self.g/((self.N0*self.N0 - self.f0*self.f0)*(sinh(2*k_z*self.Lz)/(4*k_z) - self.Lz/2)) );
                        A = (sinh(k_z*self.Lz)^2 + (self.Lz/(2*self.g))*(self.N0*self.N0 - self.f0*self.f0)*(sinh(2*k_z*self.Lz)/(2*k_z*self.Lz)-1)).^(-1/2);
                    case Normalization.omegaConstant
                        A = 1/( h0 * k_z * sqrt(1/2 + sinh(2*k_z*self.Lz)./(4*k_z*self.Lz)));
                    case Normalization.wMax
                        A = 1/sinh(k_z*self.Lz);
                    case Normalization.uMax
                        A = 1/(h0*k_z*cosh(k_z*self.Lz));
                end
                G0 = A*sinh(k_z*(self.z + self.Lz));
                F0 = A*h0*k_z*cosh(k_z*(self.z + self.Lz));
            elseif strcmp(solutionType, 'trig')
                switch self.normalization
                    case Normalization.kConstant
                        A = sqrt(self.g/((self.N0*self.N0 - self.f0*self.f0) * (self.Lz/2 - sin(2*k_z*self.Lz)/(4*k_z))));
                        A = (sin(k_z*self.Lz)^2 + (self.Lz/(2*self.g))*(self.N0*self.N0 - self.f0*self.f0)*(1-sin(2*k_z*self.Lz)/(2*k_z*self.Lz))).^(-1/2);
                    case Normalization.omegaConstant
                        A = 1/( h0 * k_z * sqrt(1/2 + sin(2*k_z*self.Lz)./(4*k_z*self.Lz)));
                    case Normalization.wMax
                        A = 1/sin(k_z*self.Lz);
                    case Normalization.uMax
                        A = 1/(h0*k_z);
                end
                G0 = A*sin(k_z*(self.z + self.Lz));
                F0 = A*h0*k_z*cos(k_z*(self.z + self.Lz));
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