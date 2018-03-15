classdef InternalModesExponentialStratification < InternalModesBase
    properties (Access = public)
        N0
        b
        rho
        N2
        rho_z
        rho_zz
        
        nInitialSearchModes = 128
        
        % analytical solutions
        shouldApproximate
        GSolutionFull
        GSolutionApprox
        FSolutionFull
        FSolutionApprox
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
            
            % A test to see which version of the solution to use             
            self.shouldApproximate = @(omega,c) (abs(bessely(b*omega/c,max(b*omega/c,(b*N0/c)*exp(-self.Lz/self.b)))/bessely(b*omega/c,(b*N0/c)*exp(-self.Lz/self.b) ))<1e-7);
            
            self.GSolutionFull = @(z,omega,c) besselj(b*(omega/c),(b*N0/c)*exp(z/self.b) )*bessely(b*(omega/c),(b*N0/c)) - bessely(b*(omega/c),(b*N0/c)*exp(z/self.b)) .* besselj(b*(omega/c),(b*N0/c));
            self.GSolutionApprox = @(z,omega,c) besselj(b*(omega/c),(b*N0/c)*exp(z/self.b) )*bessely(b*(omega/c),(b*N0/c));
            
            self.FSolutionFull = @(z,omega,c) (N0*exp(z/b)*c/2/g) .* ( (besselj(b*(omega/c)-1,(b*N0/c)*exp(z/b)) - besselj(b*(omega/c) + 1,(b*N0/c)*exp(z/b))) .*bessely(b*(omega/c),b*N0/c) - (bessely(b*(omega/c)-1,(b*N0/c)*exp(z/b))-bessely(b*(omega/c)+1,(b*N0/c)*exp(z/b))) .* besselj(b*(omega/c),b*N0/c) );
            self.FSolutionApprox = @(z,omega,c) (N0*exp(z/b)*c/2/g) .*  (besselj(b*(omega/c)-1,(b*N0/c)*exp(z/b)) - besselj(b*(omega/c) + 1,(b*N0/c)*exp(z/b))) .*bessely(b*(omega/c),b*N0/c);

            fprintf('Using the analytical form for exponential stratification N0=%.7g and b=%d\n',self.N0,self.b);
        end
                
        function [F,G,h,omega] = ModesAtWavenumber(self, k )            
            epsilon = self.f0/self.N0;
            lambda = k*self.b;
            
            x_lf = @(j,lambda) (j-1/4)*pi + lambda*pi/2;
            x_hf = @(j,lambda) lambda.*(1+0.5*(3*pi*(4*j-1)./(lambda*8*sqrt(2))).^(2/3));
            
            if lambda < 2*(1-1/4)*1e-1
                x_lowerbound = @(lambda) x_lf(1,lambda);
            else
                x_lowerbound = @(lambda) x_hf(1,lambda);
            end
            if lambda < (self.nInitialSearchModes-1/4)
                x_upperbound = @(lambda) x_lf(self.nInitialSearchModes*1.1,lambda);
            else
                x_upperbound = @(lambda) x_hf(self.nInitialSearchModes*5,lambda);
            end
            
            bounds = [x_lowerbound(lambda) x_upperbound(lambda)];
            r = FindRootsInRange(self, epsilon, lambda, bounds);
            
            while length(r) < self.nModes
                % the roots get closer together
                dr = r(end)-r(end-1);
                bounds = [bounds(2) bounds(2)+dr*self.nInitialModes];
                more_roots = FindRootsInRange(self, epsilon, lambda, bounds);
                r = [r; more_roots];
            end
            
            h = reshape((self.b*self.N0./r).^2/self.g,1,[]);
            h = h(1:self.nModes);
            
            omega = self.omegaFromK(h,k);
            
            [F,G] = NormalizedModesForOmegaAndC(self,omega,sqrt(self.g*h));
        end
        
        function r = FindRootsInRange(self, epsilon, lambda, bounds)
            x = linspace(bounds(1),bounds(2),self.nInitialSearchModes); % the choice of nInitialModes is somewhat arbitrary here.
            
            omega = @(x) sqrt( epsilon^2 * x.^2 + lambda^2 );
            if self.upperBoundary == UpperBoundary.rigidLid
                A = @(x) bessely(omega(x),x);
                B = @(x) - besselj(omega(x),x);
            elseif self.upperBoundary == UpperBoundary.freeSurface
                alpha = self.b*self.N0*self.N0/(2*self.g);
                A = @(x) bessely(omega(x),x) - (alpha./x) .* ( bessely(omega(x)-1,x) - bessely(omega(x)+1,x) );
                B = @(x) - besselj(omega(x),x) + (alpha./x) .* ( besselj(omega(x)-1,x) - besselj(omega(x)+1,x) );
            end
            f_smallnu = @(x) A(x) .* besselj(omega(x),exp(-self.Lz/self.b)*x) + B(x) .* bessely(omega(x),exp(-self.Lz/self.b)*x);
            f_bignu = @(x) (A(x) ./ bessely(omega(x),exp(-self.Lz/self.b)*x) ) .* besselj(omega(x),exp(-self.Lz/self.b)*x) + B(x);
            
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
              
            % This is the function that we use to find the eigenvalues,
            % by finding its roots.
            if omega > self.N0*exp(-self.Lz/self.b)
                eta = (sqrt(self.N0^2 - omega^2) - omega*acos(omega/self.N0))/pi;
                bounds = [0.5 self.nModes+1]; % the WKB solution should range from [3/4 nModes-1/4]
                omega_bar = omega/eta;
                N_bar = self.N0/eta;
                
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
            
            [F,G] = NormalizedModesForOmegaAndC(self,omega*ones(size(h)),sqrt(self.g*h));
            
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
        
        function [F,G] = ModeFunctionsForOmegaAndC(self,omega,c)
            if self.shouldApproximate(omega,c) == 1
                G = self.GSolutionApprox;
                F = self.FSolutionApprox;
            else
                G = self.GSolutionFull;
                F = self.FSolutionFull;
            end
        end
        
        function [F,G] = NormalizedModesForOmegaAndC(self,omega,c)
            F = zeros(length(self.z),self.nModes);
            G = zeros(length(self.z),self.nModes);
            
            for j=1:length(c)
                lowerIntegrationBound = max(5*self.b*log(omega(j)/self.N0),-self.Lz);
                [Ffunc,Gfunc] = self.ModeFunctionsForOmegaAndC(omega(j),c(j));
                switch self.normalization
                    case Normalization.uMax
                        % Doesn't the surface have the maximum value???
                        A = Ffunc(0,omega(j),c(j));
                    case Normalization.wMax
                        A = max( abs( Gfunc(linspace(lowerIntegrationBound,0,5000), omega(j), c(j) ) )  );
                    case Normalization.kConstant
                        A = sqrt(integral( @(z) (self.N0^2*exp(2*z/self.b) - self.f0^2).*Gfunc(z,omega(j),c(j)).^2,lowerIntegrationBound,0)/self.g);
                    case Normalization.omegaConstant
                        A = sqrt(integral( @(z) Ffunc(z,omega(j),c(j)).^2,lowerIntegrationBound,0)/self.Lz);
                end
                if Ffunc(0,omega(j),c(j)) < 0
                    A = -A;
                end
                F(:,j) = Ffunc(self.z,omega(j),c(j))/A;
                G(:,j) = Gfunc(self.z,omega(j),c(j))/A;
            end
        end
        
    end
    
    methods (Access = protected)
        function self = InitializeWithGrid(self, rho, zIn)
            if isempty(self.nModes) || self.nModes < 1
                self.nModes = 64;
            end
        end
        
        function self = InitializeWithFunction(self, rho, zMin, zMax, zOut)
            if isempty(self.nModes) || self.nModes < 1
                self.nModes = 64;
            end
        end
    end
    
end