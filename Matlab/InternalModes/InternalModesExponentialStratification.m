classdef InternalModesExponentialStratification < InternalModesBase
    properties (Access = public)
        N0
        b
        rho
        N2
        rho_z
        rho_zz
    end
    
    methods
        function self = InternalModesExponentialStratification(rho, z_in, z_out, latitude, varargin) 
            % rho should be a two component vector with the buoyancy at the
            % surface and the e-fold scale, e.g., [5.2e-3 1300].
            if isa(rho,'numeric') == true && length(rho) == 2
                N0 = rho(1);
                b = rho(2);
                rho0 = 1025;
                g = 9.81;
                rhoFunction = @(z) rho0*(1 + b*N0*N0/(2*g)*(1 - exp(2*z/b)));
                N2Function = @(z) N0*N0*exp(2*z/b);
            else
                error('Invalid initialization: rho must be a two-component vector with the bouyancy at the surface and the e-fold scale, e.g., [5.2e-3 1300].\n');
            end
            
            self@InternalModesBase(rhoFunction,z_in,z_out,latitude, varargin{:});
            self.N0 = N0;
            self.b = b;
            
            self.rho = rhoFunction(self.z);
            self.N2 = N2Function(self.z);
            self.rho_z = -(self.rho0*self.N0*self.N0/self.g)*exp(2*self.z/self.b);
            self.rho_zz = -(2*self.rho0*self.N0*self.N0/self.g/self.b)*exp(2*self.z/self.b);
            
            fprintf('Using the analytical form for exponential stratification N0=%.7g and b=%d\n',self.N0,self.b);
        end
                
        function [F,G,h,omega] = ModesAtWavenumber(self, k )
            error('Not yet implemented.');
            omega = self.omegaFromK(h,k);
        end
        
        function [F,G,h,k] = ModesAtFrequency(self, omega )
            if self.upperBoundary == UpperBoundary.freeSurface
                error('Not yet implemented.');
            end
            
            if omega > self.N0*exp(-self.Lz/self.b)
                eta = (sqrt(self.N0^2 - omega^2) - omega*acos(omega/self.N0))/pi;
                bounds = [0.5 self.nModes+1]; % the WKB solution should range from [3/4 nModes-1/4]
                omega_bar = omega/eta;
                N_bar = self.N0/eta;
                zT = self.b*log(omega/self.N0);
                f = @(x) besselj(omega_bar.*x,N_bar*exp(-self.Lz/self.b).*x).*bessely(omega_bar.*x,N_bar.*x)./bessely(omega_bar.*x,N_bar*exp(-self.Lz/self.b).*x) - besselj(omega_bar.*x,N_bar.*x);
                Gfunc = @(z,c) besselj(omega*self.b./c,self.N0*self.b*exp(z/self.b)./c)*bessely(omega*self.b./c,self.N0*self.b./c) - (abs(bessely(omega*self.b./c,self.N0*self.b*exp(zT/self.b)./c)/bessely(omega*self.b./c,self.N0*self.b*exp(-self.Lz/self.b)./c))>1e-7) * bessely(omega*self.b./c,self.N0*self.b*exp(z/self.b)./c) .* besselj(omega*self.b./c,self.N0*self.b./c);
                
%                 f = @(x) besselj(omega_bar.*x,N_bar*exp(-self.Lz/self.b).*x).* ( besselj(omega.*x,N_bar.*x).*cos(omega.*x*pi)-besselj(-omega.*x,N_bar.*x) ) - besselj(omega_bar.*x,N_bar.*x).*( besselj(omega.*x,N_bar*exp(-self.Lz/self.b).*x).*cos(omega.*x*pi)-besselj(-omega.*x,N_bar*exp(-self.Lz/self.b).*x) );
%                 Gfunc = @(z,c) besselj(omega*self.b./c,self.N0*self.b*exp(z/self.b)./c).* ( besselj(omega*self.b./c,self.N0*self.b./c).*cos(omega*self.b./c*pi)-besselj(-omega*self.b./c,self.N0*self.b./c) ) - besselj(omega*self.b./c,self.N0*self.b./c).*( besselj(omega*self.b./c,self.N0*self.b*exp(z/self.b)./c).*cos(omega*self.b./c*pi)-besselj(-omega*self.b./c,self.N0*self.b*exp(-z/self.b)./c) );

            else
                eta = (sqrt(self.N0^2 - omega^2) - sqrt(self.N0^2*exp(-2*self.Lz/self.b) - omega^2) - omega*acos(omega/self.N0) + omega*acos(omega/self.N0*exp(self.Lz/self.b)))/pi;
                bounds = [0.5 self.nModes+1]; % the WKB solution should range from [1 nModes]
                omega_bar = omega/eta;
                N_bar = self.N0/eta;
                f = @(x) besselj(omega_bar.*x,N_bar*exp(-self.Lz/self.b).*x).*bessely(omega_bar.*x,N_bar.*x) - besselj(omega_bar.*x,N_bar.*x).*bessely(omega_bar.*x,N_bar*exp(-self.Lz/self.b).*x);
                Gfunc = @(z,c) besselj(omega*self.b./c,self.N0*self.b*exp(z/self.b)./c) .* bessely(omega*self.b./c,self.N0*self.b./c) - besselj(omega*self.b./c,self.N0*self.b./c).*bessely(omega*self.b./c,self.N0*self.b*exp(z/self.b)./c);
            end
            
            % sqrt(gh)= b*eta/x, so h=(b*eta/x)^2/g
            f_cheb = chebfun(f,bounds,'splitting','on');
            r = roots(f_cheb);
            h = reshape((self.b*eta./r).^2/self.g,1,[]);
            h = h(1:self.nModes);
            
            F = zeros(length(self.z),self.nModes);
            G = zeros(length(self.z),self.nModes);
            
            bounds = [-self.Lz 0];
            w = chebfun(@(z) self.N0^2*exp(2*z/self.b) - self.f0^2,bounds);
            for j=1:self.nModes
                Gtmp = @(z) Gfunc(z,sqrt(self.g*h(j)));
                Gj = chebfun(Gtmp,bounds,'splitting','on');
                Fj = h(j)*diff(Gj);
                switch self.normalization
                    case Normalization.uMax
                        r = roots(diff(Fj));
                        A = max( abs( Fj(r) ) );
                    case Normalization.wMax
                        r = roots(Fj);
                        A = max( abs( Gj(r) ) );
                    case Normalization.kConstant           
                        A = sqrt(sum(w.*Gj.^2)/self.g);
                    case Normalization.omegaConstant
                        A = sqrt(sum(Fj.^2)/self.Lz);
                end
                if Fj(0) < 0
                    A = -A;
                end
                F(:,j) = Fj(self.z)/A;
                G(:,j) = Gj(self.z)/A;
            end
            
            k = self.kFromOmega(h,omega);
        end
        
        function [psi] = SurfaceModesAtWavenumber(self, k)
            % size(psi) = [size(k); length(z)]
            error('Not yet implemented. See LaCasce 2012.');
        end
        
        function [psi] = BottomModesAtWavenumber(self, k)
            % size(psi) = [size(k); length(z)]
            error('Not yet implemented. See LaCasce 2012.');
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
    
end