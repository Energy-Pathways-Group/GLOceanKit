classdef InternalModesExponentialStratification < InternalModesBase
    properties (Access = public)
        N0
        b
        rho
        N2
        rho_z
        rho_zz
        
        om
        s
        shouldTruncate
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
            
            % omega/(omega*exp(-self.Lz/self.b)
            % omega/(self.N0*exp(-self.Lz/self.b)
            self.om = @(omega,c) omega*self.b./c;
            self.s = @(z,c) self.N0*self.b*exp(z/self.b)./c;
            self.shouldTruncate = @(omega,c,zT) (abs(bessely(self.om(omega,c),self.s(zT,c))/bessely(self.om(omega,c),self.s(-self.Lz,c) ))>1e-7);
            
            fprintf('Using the analytical form for exponential stratification N0=%.7g and b=%d\n',self.N0,self.b);
        end
                
        function [F,G,h,omega] = ModesAtWavenumber(self, k )            
            epsilon = self.f0/self.N0;
            lambda = k*self.b;
            
            x_lf = @(j,lambda) (j-1/4)*pi + lambda*pi/2;
            x_hf = @(j,lambda) lambda.*(1+0.5*(3*pi*(4*j-1)./(lambda*8*sqrt(2))).^(2/3));
            
            self.nInitialModes = 128;
            x_lowerbound = @(lambda) (lambda < 2*(1-1/4)*1e-1).*x_lf(1,lambda) + (lambda >= 2*(1-1/4)*1e-1).*x_hf(1,lambda);
            x_upperbound = @(lambda) (lambda < (nInitialModes-1/4)).*x_lf(nInitialModes*1.1,lambda) + (lambda >= (nInitialModes-1/4)).*x_hf(nInitialModes*5,lambda);
            
            bounds = [x_lowerbound(lambda) x_upperbound(lambda)];
            r = FindRootsInRange(self, epsilon, lambda, bounds);
            
            while length(r) < self.nModes
                % the roots get closer together
                dr = r(end)-r(end-1);
                bounds = [bounds(2) bounds(2)+dr*self.nInitialModes];
                more_roots = FindRootsInRange(self, epsilon, lambda, bounds);
                r = [r; more_roots];
            end
            
            h = reshape((self.b*eta./r).^2/self.g,1,[]);
            h = h(1:self.nModes);
            
            omega = self.omegaFromK(h,k);
        end
        
        function r = FindRootsInRange(self, epsilon, lambda, bounds)
            x = linspace(bounds(1),bounds(2),self.nInitialModes); % the choice of nInitialModes is somewhat arbitrary here.
            
            omega = @(x) sqrt( epsilon^2 * x.^2 + lambda^2 );
            f_smallnu = @(x) besselj(omega(x),exp(-self.Lz/self.b)*x).*bessely(omega(x),x) - besselj(omega(x),x).*bessely(omega(x),exp(-self.Lz/self.b)*x);
            f_bignu = @(x) besselj(omega(x),exp(-self.Lz/self.b)*x).*bessely(omega(x),x)./bessely(omega(x),exp(-self.Lz/self.b)*x) - besselj(omega(x),x);
            
            % The function omega(x)./(exp(-D/b)*x) will monotonically decay with x.
            % We want to find where it first drops below 5, and use the appropriate
            % form of the equation.
            xcutoffIndex = find( omega(x)./(exp(-self.Lz/self.b)*x) < 5,1,'first');
            
            if xcutoffIndex == 1
                % nu is small for all values of x
                f = f_smallnu;
                f_cheb = chebfun(f,bounds,'splitting','on');
                r = roots(f_cheb);
            elseif isempty(xcutoffIndex) == 1
                % nu is large for all values of x
                f = f_bignu;
                f_cheb = chebfun(f,bounds,'splitting','on');
                r = roots(f_cheb);
            else
                fprintf('big nu/small nu\n');
                % we need to subdivide the interval into small and large.
                f = f_bignu;
                f_cheb = chebfun(f,[bounds(1) x(xcutoffIndex)] ,'splitting','on');
                r = roots(f_cheb);
                f = f_smallnu;
                f_cheb = chebfun(f,[x(xcutoffIndex) bounds(2)] ,'splitting','on');
                r = [r; roots(f_cheb)];
            end
        end
        
        function [F,G,h,k] = ModesAtFrequency(self, omega )
            if self.upperBoundary == UpperBoundary.freeSurface
                error('Not yet implemented.');
            end
            
            zT = max(self.b*log(omega/self.N0),-self.Lz);
            if omega > self.N0*exp(-self.Lz/self.b)
                eta = (sqrt(self.N0^2 - omega^2) - omega*acos(omega/self.N0))/pi;
                bounds = [0.5 self.nModes+1]; % the WKB solution should range from [3/4 nModes-1/4]
                omega_bar = omega/eta;
                N_bar = self.N0/eta;
                
                
                % This is the function that we use to find the eigenvalues,
                % by finding its roots.
                f = @(x) besselj(omega_bar.*x,N_bar*exp(-self.Lz/self.b).*x).*bessely(omega_bar.*x,N_bar.*x)./bessely(omega_bar.*x,N_bar*exp(-self.Lz/self.b).*x) - besselj(omega_bar.*x,N_bar.*x);
                
            else
                eta = (sqrt(self.N0^2 - omega^2) - sqrt(self.N0^2*exp(-2*self.Lz/self.b) - omega^2) - omega*acos(omega/self.N0) + omega*acos(omega/self.N0*exp(self.Lz/self.b)))/pi;
                bounds = [0.5 self.nModes+1]; % the WKB solution should range from [1 nModes]
                omega_bar = omega/eta;
                N_bar = self.N0/eta;
                f = @(x) besselj(omega_bar.*x,N_bar*exp(-self.Lz/self.b).*x).*bessely(omega_bar.*x,N_bar.*x) - besselj(omega_bar.*x,N_bar.*x).*bessely(omega_bar.*x,N_bar*exp(-self.Lz/self.b).*x);
            end
            
            % sqrt(gh)= b*eta/x, so h=(b*eta/x)^2/g
            f_cheb = chebfun(f,bounds,'splitting','on');
            r = roots(f_cheb);
            h = reshape((self.b*eta./r).^2/self.g,1,[]);
            h = h(1:self.nModes);
            
            F = zeros(length(self.z),self.nModes);
            G = zeros(length(self.z),self.nModes);
            
            lowerIntegrationBound = max(5*zT,-self.Lz);
            for j=1:self.nModes
                c = sqrt(self.g*h(j));
                [Ffunc,Gfunc] = self.ModeFunctionsForOmegaAndC(omega,c,zT);
                switch self.normalization
                    case Normalization.uMax
                        A = max( abs( Ffunc(linspace(lowerIntegrationBound,0,5000),c ) ) );
                    case Normalization.wMax
                        A = max( abs( Gfunc(linspace(lowerIntegrationBound,0,5000),c ) )  );
                    case Normalization.kConstant           
                        A = sqrt(integral( @(z) (self.N0^2*exp(2*z/self.b) - self.f0^2).*Gfunc(z,c).^2,lowerIntegrationBound,0)/self.g);
                    case Normalization.omegaConstant
                        A = sqrt(integral( @(z) Ffunc(z,c).^2,lowerIntegrationBound,0)/self.Lz);
                end
                if Ffunc(0,c) < 0
                    A = -A;
                end
                F(:,j) = Ffunc(self.z,c)/A;
                G(:,j) = Gfunc(self.z,c)/A;
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
        
        function [F,G] = ModeFunctionsForOmegaAndC(self,omega,c,zT)
            if nargin == 2
                zT = -self.Lz;
            end
            omega_bar = @(c) self.om(omega,c);
            s_bar = self.s;
            %if self.shouldTruncate(omega,c,zT) == 0
            if (omega/self.N0)*exp(self.Lz/self.b) > 5
                fprintf('Using truncated solution.\n');
                G = @(z,c) besselj( omega_bar(c),s_bar(z,c) )*bessely(omega_bar(c),s_bar(0,c));
                F = @(z,c) (self.N0*exp(z/self.b)*c/2/self.g) .* ( (besselj(omega_bar(c)-1,s_bar(z,c)) - besselj(omega_bar(c) + 1,s_bar(z,c))) .*bessely(omega_bar(c),s_bar(0,c)) );
            else
                fprintf('Using full solution.\n');
                G = @(z,c) besselj( omega_bar(c),s_bar(z,c) )*bessely(omega_bar(c),s_bar(0,c)) - bessely(omega_bar(c),s_bar(z,c)) .* besselj(omega_bar(c),s_bar(0,c));
                F = @(z,c) (self.N0*exp(z/self.b)*c/2/self.g) .* ( (besselj(omega_bar(c)-1,s_bar(z,c)) - besselj(omega_bar(c) + 1,s_bar(z,c))) .*bessely(omega_bar(c),s_bar(0,c)) - (bessely(omega_bar(c)-1,s_bar(z,c))-bessely(omega_bar(c)+1,s_bar(z,c))) .* besselj(omega_bar(c),s_bar(0,c)) );
            end
        end
        
        function G = GModeForOmegaAndC(self,omega,c)
            % omega/(omega*exp(-self.Lz/self.b)
            % omega/(self.N0*exp(-self.Lz/self.b)
%             self.om = @(omega,c) omega*self.b./c;
%             self.s = @(z,c) self.N0*self.b*exp(z/self.b)./c;
% zT = max(self.b*log(omega/self.N0),-self.Lz);
            nu = omega*self.b/c;
            x_max = self.N0*self.b/c;
            x_min = x_max*exp(-self.Lz/self.b);
            x = max(x_min, nu);
            if abs(bessely(nu,x)/bessely(nu,x_min )) > 1e-7
                G = @(z,c) besselj(nu,x_max*exp(z/self.b) )*bessely(nu,x_max);
            else
                G = @(z,c) besselj(nu,x_max*exp(z/self.b) )*bessely(nu,x_max) - bessely(nu,x_max*exp(z/self.b)) .* besselj(nu,x_max);
            end
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