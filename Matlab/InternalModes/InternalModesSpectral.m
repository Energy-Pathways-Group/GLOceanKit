classdef InternalModesSpectral < InternalModesBase
    % InternalModesSpectral uses Chebyshev polynomials on a z-grid to
    % compute the internal wave modes. See InternalModesBase for basic
    % usage information.
    %
    % This class takes the name/value pair 'nEVP' (default 513) in order to
    % set the default resolution of the polynomials used to solve the
    % eigenvalue problem (EVP), e.g.,
    %       modes = InternalModes(rho,zDomain,zOut,latitude, 'nEVP', 128);
    % Generally speaking, you will get half or more good quality modes (so
    % nEVP=513 should give you 250 or so quality modes, if the density
    % structure isn't too weird). Note that the time to solve the EVP
    % scales as the nEVP^3, so going higher resolution comes at great cost.
    % For best speed results, nEVP should be set to 2^n+1 to fully utilize
    % the FFT.
    %
    % All solutions are spectrally projected to the request ouput grid, and
    % therefore will remain high quality. The fastest output can be
    % achieved by using an Gauss-Lobatto grid spanning z.
    % 
    % The modes are normalized by integrating the Chebyshev polynomials on
    % the Lobatto grid using the exact integrals.
    %
    % A word on notation:
    %
    % Dimensions -- like z -- have grids, which may differ. zIn is the
    % dimension z, on some grid given by the user. zLobatto is the
    % diension z, on a Lobatto grid. zCheb would be used for the Chebyshev
    % polynomial representation of a function defined on z. Note that the
    % base class property 'z' is really zOut.
    %
    % Functions -- if they're analytical, they're defined on a dimension,
    % not a grid. If they're gridded, then the grid must also be specified.
    % Thus, the function rho(z) when analytically defined will be given as
    % either just rho (when obvious) or rho_z when necessary. If it's
    % gridded, then expect rho_zIn or rho_zLobatto.
    %
    % Transformations -- a transformation often a matrix (or generally just
    % a function) that transformation a function from one basis to another.
    % A common scenario here is that we need to go from zLobatto to zCheb.
    % In our cases we're usually going from one gridded dimension to
    % another. Inverse Chebyshev transformations are denote by T, or more
    % specifically T_zCheb_zLobatto.
    %
    % The underscore in the case of rho_z is certainly recognized as a
    % derivative (because of LaTex notation), so it is awkward with the
    % notion used here, where the underscore is used to separate the
    % function/transformation name from the grid. Thus, we leave rho_z and
    % rho_zz as the only two exceptions to the aforementioned notation
    % because they are the only public facing interface.
    %
    % The 'x' dimension used in the properties for this class are denoted
    % as such because the different subclasses uses these properties to
    % store values on different grids. In other words, this class uses 'x'
    % for the z dimension (depth), but the density subclass uses 'x' for
    % its density stretched coordinated.
    %
    %   See also INTERNALMODES, INTERNALMODESBASE,
    %   INTERNALMODESDENSITYSPECTRAL, INTERNALMODESWKBSPECTRAL, and
    %   INTERNALMODESFINITEDIFFERENCE.
    %
    %
    %   Jeffrey J. Early
    %   jeffrey@jeffreyearly.com
    %
    %   March 14th, 2017        Version 1.0
    
    properties (Access = public)

    end
    
    properties (Dependent)
        rho % Density on the z grid.
        rho_z % First derivative of density on the z grid.
        rho_zz % Second derivative of density on the z grid.
        N2 % Buoyancy frequency on the z grid, $N^2 = -\frac{g}{\rho(0)} \frac{\partial \rho}{\partial z}$.
    end
    

    properties %(Access = private)
        rho_function        % function handle to return rho at z
        N2_function         % function handle to return N2 at z
                
        % These properties are initialized with SetupEigenvalueProblem()
        % Most subclasses *will* override these initializations. The 'x' refers to the stretched coordinate being used.
        % This class uses x=z (depth), although they may have different numbers of points.
        nEVP = 0           % number of points in the eigenvalue problem
        xLobatto           % stretched coordinate Lobatto grid nEVP points. z for this class, density or wkb for others.
        xDomain            % limits of the stretched coordinate [xMin xMax]
        z_xLobatto         % The value of z, at the xLobatto points
        xOut               % desired locations of the output in x-coordinate (deduced from z_out)
        N2_xLobatto        % N2 on the z2Lobatto grid
        Diff1_xCheb        % single derivative in spectral space, *function handle*
        T_xLobatto, Tx_xLobatto, Txx_xLobatto        % Chebyshev polys (and derivs) on the zLobatto
        T_xCheb_zOut
        Int_xCheb           % Vector that multiplies Cheb coeffs, then sum for integral
    end
    
    properties (Dependent)
        xMin
        xMax
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = InternalModesSpectral(rho, zIn, zOut, latitude,varargin)
            self@InternalModesBase(rho,zIn,zOut,latitude,varargin{:});
            self.SetupEigenvalueProblem();            
        end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computed (dependent) properties
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function value = get.rho(self)
            value = self.rho_function(self.z);
        end
                
        function value = get.rho_z(self)
            rho_z_function = diff(self.rho_function);
            value = rho_z_function(self.z);
        end
        
        function value = get.rho_zz(self)
            rho_zz_function = diff(diff(self.rho_function));
            value = rho_zz_function(self.z);
        end
        
        function value = get.N2(self)
            value = self.N2_function(self.z);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computation of the modes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [A,B] = EigenmatricesForWavenumber(self, k )
            % The eigenvalue equation is,
            % G_{zz} - K^2 G = \frac{f_0^2 -N^2}{gh_j}G
            % A = \frac{g}{f_0^2 -N^2} \left( \partial_{zz} - K^2*I \right)
            % B = I
            T = self.T_xLobatto;
            Tz = self.Tx_xLobatto;
            Tzz = self.Txx_xLobatto;
            n = self.nEVP;
            
            A = (Tzz - k*k*eye(n)*T);
            B = diag((self.f0*self.f0-self.N2_xLobatto)/self.g)*T;
            
            switch self.lowerBoundary
                case LowerBoundary.rigidLid
                    A(n,:) = T(n,:);
                    B(n,:) = 0;
                case LowerBoundary.none
                otherwise
                    error('Unknown boundary condition');
            end
            
            switch self.upperBoundary
                case UpperBoundary.freeSurface
                    % G_z = \frac{1}{h_j} G at the surface
                    A(1,:) = Tz(1,:);
                    B(1,:) = T(1,:);
                case UpperBoundary.rigidLid
                    A(1,:) = T(1,:);
                    B(1,:) = 0;
                case UpperBoundary.none
                otherwise
                    error('Unknown boundary condition');
            end
        end
        
        function [A,B] = EigenmatricesForFrequency(self, omega )
            T = self.T_xLobatto;
            Tz = self.Tx_xLobatto;
            Tzz = self.Txx_xLobatto;
            n = self.nEVP;
            
            A = Tzz;
            B = diag((omega*omega-self.N2_xLobatto)/self.g)*T;
            
            switch self.lowerBoundary
                case LowerBoundary.rigidLid
                    A(n,:) = T(n,:);
                    B(n,:) = 0;
                case LowerBoundary.none
                otherwise
                    error('Unknown boundary condition');
            end
            
            switch self.upperBoundary
                case UpperBoundary.freeSurface
                    % G_z = \frac{1}{h_j} G at the surface
                    A(1,:) = Tz(1,:);
                    B(1,:) = T(1,:);
                case UpperBoundary.rigidLid
                    A(1,:) = T(1,:);
                    B(1,:) = 0;
                case UpperBoundary.none
                otherwise
                    error('Unknown boundary condition');
            end
        end
        
        function [F,G,h,omega,F2,N2G2] = ModesAtWavenumber(self, k )
            self.gridFrequency = 0;
            
            [A,B] = self.EigenmatricesForWavenumber(k);
            
            if nargout == 6
                [F,G,h,F2,N2G2] = self.ModesFromGEPSpectral(A,B);
            else
                [F,G,h] = self.ModesFromGEPSpectral(A,B); 
            end
            omega = self.omegaFromK(h,k);
        end
        
        function [F,G,h,k,F2,N2G2] = ModesAtFrequency(self, omega )
            self.gridFrequency = omega;
            
            [A,B] = self.EigenmatricesForFrequency(omega);
            
            if nargout == 6
                [F,G,h,F2,N2G2] = self.ModesFromGEPSpectral(A,B);
            else
                [F,G,h] = self.ModesFromGEPSpectral(A,B);
            end
            k = self.kFromOmega(h,omega);
        end 
        
        function psi = SurfaceModesAtWavenumber(self, k) 
            psi = self.BoundaryModesAtWavenumber(k,1);
        end
        
        function psi = BottomModesAtWavenumber(self, k) 
            psi = self.BoundaryModesAtWavenumber(k,0);
        end
        
        function psi = BoundaryModesAtWavenumber(self, k, isSurface)            
            % Estimate the grid resolution necessary to resolve the
            % smallest mode.
            sizeK = size(k);
            if length(sizeK) == 2 && sizeK(2) == 1
                sizeK(2) = [];
            end
            K = k(:);
            maxK = max(K(K>0));
            maxN = sqrt(max(self.N2(1),self.N2(end)));
            minDz = self.f0/(maxK*maxN)/10;
            n = max(128,2^ceil(log2((pi/2)*sqrt(self.Lz/minDz)+1))+1);
%             fprintf('Using %d grid points for the SQG mode\n',n);
            
            % Now create this new grid
            yLobatto = (self.Lz/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + min(self.zLobatto);
            T_zCheb_yLobatto = InternalModesSpectral.ChebyshevTransformForGrid(self.zLobatto, yLobatto);
            T_yCheb_zOut = InternalModesSpectral.ChebyshevTransformForGrid(yLobatto, self.z);
            
            [T,Ty,Tyy] = InternalModesSpectral.ChebyshevPolynomialsOnGrid( yLobatto, length(yLobatto) );
            N2_yLobatto = T_zCheb_yLobatto(self.N2_zCheb);
            N2_z_yLobatto = T_zCheb_yLobatto(self.Diff1_zCheb(self.N2_zCheb));
            
            A = N2_yLobatto .* Tyy - N2_z_yLobatto.*Ty;
            B = - (1/(self.f0*self.f0))* (N2_yLobatto.*N2_yLobatto) .*T;
            
            b = zeros(size(yLobatto));
            if isSurface == 1
                b(1) = 1;
            else
                b(end) = 1;
            end
            
            psi = zeros(length(k),length(self.z));
            for ii = 1:length(k)
                M = A + k(ii)*k(ii)*B;
                M(1,:) = self.f0*Ty(1,:);
                M(end,:) = self.f0*Ty(end,:);
                
                psi_cheb = M\b;
                psi(ii,:) = T_yCheb_zOut(psi_cheb);
            end
            
            sizeK(end+1) = length(self.z);
            psi = reshape(psi,sizeK);
        end
        
        function z_g = GaussQuadraturePointsForModesAtWavenumber(self,nPoints,k)
            % Now we just need to find the roots of the n+1 mode.
            % For constant stratification this should give back the
            % standard Fourier modes, i.e., an evenly spaced grid.
            %
            % Note that if the boundary conditions are such that G(0)=0 and
            % G(-D)=0, then those two points do not encode any information.
            % As such, only the first (nPoints-2) modes will encode any
            % useful information. So we'd expect cond(G(:,1:(nPoints-2))))
            % to be good (low), but not the next.
            if 2*nPoints < self.nEVP
               [A,B] = self.EigenmatricesForWavenumber(k);
               if ( any(any(isnan(A))) || any(any(isnan(B))) )
                   error('EVP setup fail. Found at least one nan in matrices A and B.\n');
               end
               [V,D] = eig( A, B );
               
               hFromLambda = @(lambda) 1.0 ./ lambda;
               [h, permutation] = sort(real(hFromLambda(diag(D))),'descend');
               G_cheb=V(:,permutation);
               maxModes = ceil(find(h>0,1,'last')/2);
               
               if maxModes < (nPoints+1)
                   error('We tried, but you are gonna need more points.');
               end
               
               if 1 == 0
                   F = self.Diff1_xCheb(G_cheb(:,nPoints-1));
                   roots = FindRootsFromChebyshevVector(F(1:end-1), self.zLobatto);
                   z_g = cat(1,min(self.zLobatto),reshape(roots,[],1),max(self.zLobatto));
               else
                   if self.upperBoundary == UpperBoundary.rigidLid
                       % n-th mode has n+1 zeros (including boundaries)
                       roots = FindRootsFromChebyshevVector(G_cheb(:,nPoints-1), self.zLobatto);
                   elseif self.upperBoundary == UpperBoundary.freeSurface
                       % n-th mode has n zeros (including zero at lower
                       % boundary, and not zero at upper)
                       a = InternalModesSpectral.ValueOfFunctionAtPointOnGrid(max(self.zLobatto),self.zLobatto,G_cheb(:,nPoints-0));
                       b = InternalModesSpectral.ValueOfFunctionAtPointOnGrid(max(self.zLobatto),self.zLobatto,G_cheb(:,nPoints-1));
                       q = G_cheb(:,nPoints-0) - (a/b)*G_cheb(:,nPoints-1);
%                        t1 = InternalModesSpectral.ValueOfFunctionAtPointOnGrid(max(self.zLobatto),self.zLobatto,q);
%                        t2 = InternalModesSpectral.ValueOfFunctionAtPointOnGrid(min(self.zLobatto),self.zLobatto,q);
%                        q = G_cheb(:,nPoints-0);
                       roots = FindRootsFromChebyshevVector(q, self.zLobatto);
                   end
                   z_g = reshape(roots,[],1);
               end
               
               z_g(z_g<min(self.zLobatto)) = min(self.zLobatto);
               z_g(z_g>max(self.zLobatto)) = max(self.zLobatto);
               z_g = unique(z_g,'stable');
            else
                error('need more points');
            end
        end
        
        function value = get.xMin(self)
            value = self.xDomain(1);
        end
        
        function value = get.xMax(self)
            value = self.xDomain(2);
        end
    end
    
    methods (Access = protected)
        
        function self = SetupEigenvalueProblem(self)
            % Called during initialization after InitializeZLobattoProperties().
            % Subclasses will override this function.
            n = self.nEVP;
            self.xLobatto = (self.Lz/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + self.zMin;
            self.xLobatto(1) = self.zMax;
            self.xDomain = self.zDomain;
            
            self.N2_xLobatto = self.N2_function(self.xLobatto);
            self.Diff1_xCheb = @(v) (2/self.Lz)*InternalModesSpectral.DifferentiateChebyshevVector( v );
            
            [self.T_xLobatto,self.Tx_xLobatto,self.Txx_xLobatto] = InternalModesSpectral.ChebyshevPolynomialsOnGrid( self.xLobatto, length(self.xLobatto) );
            self.T_xCheb_zOut = InternalModesSpectral.ChebyshevTransformForGrid(self.xLobatto, self.z);
            
            % We use that \int_{-1}^1 T_n(x) dx = \frac{(-1)^n + 1}{1-n^2}
            % for all n, except n=1, where the integral is zero.
            np = (0:(n-1))';
            self.Int_xCheb = -(1+(-1).^np)./(np.*np-1);
            self.Int_xCheb(2) = 0;
            self.Int_xCheb = self.Lz/2*self.Int_xCheb;
            
            if self.shouldShowDiagnostics == 1
                fprintf(' The eigenvalue problem will be solved with %d points.\n', length(self.xLobatto));
            end
        end
        
        function self = InitializeWithGrid(self, rho, zIn)
            self.validateInitialModeAndEVPSettings();
            
            K = 6; % cubic spline
            if self.requiresMonotonicDensity == 1
                if any(diff(rho)./diff(zIn) > 0)
                    rho_interpolant = SmoothingSpline(zIn,rho,NormalDistribution(1),'K',K,'knot_dof',1,'lambda',Lambda.optimalExpected,'constraints',struct('global',ShapeConstraint.monotonicDecreasing));
                    rho_interpolant.minimize( @(spline) spline.expectedMeanSquareErrorFromCV );
                    if self.shouldShowDiagnostics == 1
                        fprintf('Creating a %d-order monotonic smoothing spline from the %d points.\n', K, length(rho));
                    end
                else
                    z_knot = InterpolatingSpline.KnotPointsForPoints(zIn,K,1);
                    rho_interpolant = ConstrainedSpline(zIn,rho,K,z_knot,NormalDistribution(1),struct('global',ShapeConstraint.monotonicDecreasing));
                    if self.shouldShowDiagnostics == 1
                        fprintf('Creating a %d-order monotonic spline from the %d points.\n', K, length(rho));
                    end
                end
            else
                rho_interpolant = InterpolatingSpline(zIn,rho,'K',K);
                if self.shouldShowDiagnostics == 1
                    fprintf('Creating a %d-order spline from the %d points.\n', K, length(rho));
                end
            end
            
            self.rho_function = rho_interpolant;
            self.N2_function = -(self.g/self.rho0)*diff(self.rho_function);
        end
        
        function self = InitializeWithFunction(self, rho, zMin, zMax)
            self.validateInitialModeAndEVPSettings();
            
            if ~exist('chebfun','class')
               error('The package chebfun is required when initializing with a function.')
            end

            self.rho_function = chebfun(rho,[zMin zMax]);
            self.N2_function = -(self.g/self.rho0)*diff(self.rho_function);
            
            if self.requiresMonotonicDensity == 1
                % Taken from inv as part of chebfun.
                f = self.rho_function;
                fp = diff(f);
                tPoints = roots(fp);
                if ( ~isempty(tPoints) )
                    endtest = zeros(length(tPoints), 1);
                    for k = 1:length(tPoints)
                        endtest(k) = min(abs(tPoints(k) - f.domain));
                    end
                    if ( any(endtest > 100*abs(feval(f, tPoints))*tol) )
                        error('Density must be monotonic to use WKB. The function you provided is not monotonic.');
                    end
                end
            end
            
            if self.shouldShowDiagnostics == 1
                fprintf('Projected the function onto %d Chebyshev polynomials\n', length(self.rho_function));
            end
        end      
                   
        function self = validateInitialModeAndEVPSettings(self)
            % The user requested that the eigenvalue problem be solved on a
            % grid of particular length
            if self.nEVP > 0
                if self.nModes > self.nEVP
                    self.nEVP = self.nModes;
                end
            else
                self.nEVP = 513; % 2^n + 1 for a fast Chebyshev transform
            end            
        end
                        
        % This function is an intermediary used by ModesAtFrequency and
        % ModesAtWavenumber to establish the various norm functions.
        function [F,G,h,F2,N2G2] = ModesFromGEPSpectral(self,A,B)
            hFromLambda = @(lambda) 1.0 ./ lambda;
            GOutFromGCheb = @(G_cheb,h) self.T_xCheb_zOut(G_cheb);
            FOutFromGCheb = @(G_cheb,h) h * self.T_xCheb_zOut(self.Diff1_xCheb(G_cheb));
            GFromGCheb = @(G_cheb,h) InternalModesSpectral.ifct(G_cheb);
            FFromGCheb = @(G_cheb,h) h * InternalModesSpectral.ifct( self.Diff1_xCheb(G_cheb) );
            GNorm = @(Gj) abs(Gj(1)*Gj(1) + sum(self.Int_xCheb .*InternalModesSpectral.fct((1/self.g) * (self.N2_xLobatto - self.f0*self.f0) .* Gj .^ 2)));
            FNorm = @(Fj) abs(sum(self.Int_xCheb .*InternalModesSpectral.fct((1/self.Lz) * Fj.^ 2)));
            if nargout == 5
                [F,G,h,F2,N2G2] = ModesFromGEP(self,A,B,hFromLambda,GFromGCheb,FFromGCheb,GNorm,FNorm,GOutFromGCheb,FOutFromGCheb);
            else
                [F,G,h] = ModesFromGEP(self,A,B,hFromLambda,GFromGCheb,FFromGCheb,GNorm,FNorm,GOutFromGCheb,FOutFromGCheb);
            end
        end
        
        % Take matrices A and B from the generalized eigenvalue problem
        % (GEP) and returns F,G,h. The last seven arguments are all
        % function handles that do as they say.
        function [F,G,h,F2,N2G2] = ModesFromGEP(self,A,B,hFromLambda,GFromGCheb, FFromGCheb, GNorm,FNorm, GOutFromGCheb,FOutFromGCheb)
            if ( any(any(isnan(A))) || any(any(isnan(B))) )
                error('EVP setup fail. Found at least one nan in matrices A and B.\n');
            end
            [V,D] = eig( A, B );
            
            [h, permutation] = sort(real(hFromLambda(diag(D))),'descend');
            G_cheb=V(:,permutation);
            
            if isempty(self.nModes)
                maxModes = ceil(find(h>0,1,'last')/2); % Have to do ceil, not floor, or we lose the barotropic mode.
                if maxModes == 0
                    fprintf('No usable modes found! Try with higher resolution.\n');
                    return;
                end
            else
                maxModes = self.nModes;
            end
            
            
            F = zeros(length(self.z),maxModes);
            G = zeros(length(self.z),maxModes);
            h = reshape(h(1:maxModes),1,[]);
            
            % only used if nargout == 5
            N2G2 = zeros(1,maxModes);
            F2 = zeros(1,maxModes);
            
            % This still need to be optimized to *not* do the transforms
            % twice, when the EVP grid is the same as the output grid.
            [maxIndexZ] = find(self.N2_xLobatto-self.gridFrequency*self.gridFrequency>0,1,'first');
            if maxIndexZ > 1 % grab a point just above the turning point, which should have the right sign.
                maxIndexZ = maxIndexZ-1;
            elseif isempty(maxIndexZ)
                maxIndexZ = 1;
            end
            for j=1:maxModes
                Fj = FFromGCheb(G_cheb(:,j),h(j));
                Gj = GFromGCheb(G_cheb(:,j),h(j));
                switch self.normalization
                    case Normalization.uMax
                        A = max( abs( Fj ));
                    case Normalization.wMax
                        A = max( abs( Gj ) );
                    case Normalization.kConstant
                        A = sqrt(GNorm( Gj ));
                    case Normalization.omegaConstant
                        A = sqrt(FNorm( Fj ));
                end
                if Fj(maxIndexZ) < 0
                    A = -A;
                end
                
                G(:,j) = GOutFromGCheb(G_cheb(:,j),h(j))/A;
                F(:,j) = FOutFromGCheb(G_cheb(:,j),h(j))/A;
                
                if nargout == 5
                    F2(j) = self.Lz*FNorm( Fj/A );
                    N2G2(j) = self.g*(GNorm( Gj/A )-Gj(1)*Gj(1)) + self.f0*self.f0*self.Lz*FNorm( Gj/A ); % this is being clever, but should give \int N2*G2 dz
                end   
            end
            

        end
        
        function zTPs = FindTurningPointsAtFrequency(self, omega)
            f_cheb = self.N2_zCheb;
            f_cheb(1) = f_cheb(1) - omega*omega;
            zTPs= InternalModesSpectral.FindRootsFromChebyshevVector(f_cheb,self.zLobatto);
        end
    end
    
    methods (Static)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Methods to find the turning points (used for normalization)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [zBoundariesAndTPs, thesign, boundaryIndices] = FindTurningPointBoundariesAtFrequency(N2, z, omega)
            % This function returns not just the turning points, but also
            % the top and bottom boundary locations in z. The boundary
            % indices are the index to the point just *above* the turning
            % point.
            N2Omega2 = N2 - omega*omega;
            a = N2Omega2; a(a>=0) = 1; a(a<0) = 0;
            turningIndices = find(diff(a)~=0);
            nTP = length(turningIndices);
            zTP = zeros(nTP,1);
            for i=1:nTP
                fun = @(depth) interp1(z,N2Omega2,depth,'spline');
                z0 = [z(turningIndices(i)+1) z(turningIndices(i))];
                zTP(i) = fzero(fun,z0);
            end
            boundaryIndices = [1; turningIndices; length(z)];
            zBoundariesAndTPs = [z(1); zTP; z(end)];
            
            % what's the sign on the EVP in these regions? In this case,
            % positive indicates oscillatory, negative exponential decay
            midZ = zBoundariesAndTPs(1:end-1) + diff(zBoundariesAndTPs)/2;
            thesign = sign( interp1(z,N2Omega2,midZ,'spline') );
            % remove any boundaries of zero length
            for index = reshape(find( thesign == 0),1,[])
                thesign(index) = [];
                zBoundariesAndTPs(index) = [];
                boundaryIndices(index) = [];
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Convert to a stretched grid
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [z_xLobatto, xOut] = StretchedGridFromCoordinate( x, xLobatto, zLobatto, zOut)
            % x is a function handle, that maps from z.
            % xLobatto is the grid on x
            % zLobatto is the Lobatto grid on z
            % zOut is the output grid
            
            maxZ = max(zLobatto);
            minZ = min(zLobatto);
%             z_xLobatto = interp1(x(zLobatto), zLobatto, xLobatto, 'spline','extrap');
            z_xLobatto = InternalModesSpectral.fInverseBisection(x,xLobatto,min(zLobatto),max(zLobatto),1e-12);
            z_xLobatto(z_xLobatto>maxZ) = maxZ;
            z_xLobatto(z_xLobatto<minZ) = minZ;
            
            s_z_xLobatto = x(z_xLobatto); % this should give us xLobatto back, if our approximation is good.
            relativeError = max( abs(s_z_xLobatto - xLobatto))/max(abs(xLobatto));
            nloops = 0;
            maxLoops = 100;
            while ( nloops < maxLoops && relativeError > 1e-10)
                z_xLobatto = interp1(s_z_xLobatto, z_xLobatto, xLobatto, 'spline','extrap');
                z_xLobatto(z_xLobatto>maxZ) = maxZ;
                z_xLobatto(z_xLobatto<minZ) = minZ;
                s_z_xLobatto = x(z_xLobatto);
            
                relativeError = max( abs(s_z_xLobatto - xLobatto))/max(abs(xLobatto));
                nloops = nloops + 1;
            end
            if nloops == maxLoops
                warning('Stretched coordinate reached maximum number of loops (%d) with relative error of %g\n', maxLoops, relativeError);
            end
                        
            xOut = x(zOut);
            xOut(xOut>max(xLobatto)) = max(xLobatto);
            xOut(xOut<min(xLobatto)) = min(xLobatto);
        end
        
        function [flag, dTotalVariation, rho_zCheb, rho_zLobatto, rhoz_zCheb, rhoz_zLobatto] = CheckIfReasonablyMonotonic(zLobatto, rho_zCheb, rho_zLobatto, rhoz_zCheb, rhoz_zLobatto)
            % We want to know if the density function is decreasing as z
            % increases. If it's not, are the discrepencies small enough
            % that we can just force them?
            %
            % flag is 0 if the function is not reasonably monotonic, 1 if
            % it is, and 2 if we were able to to coerce it to be, without
            % too much change
            
            flag = 0;
            dTotalVariation = 0;
            if any(rhoz_zLobatto > 0)
                % record the density at the bottom, and the total variation
                % in density
                rho_top = min(rho_zLobatto(1),rho_zLobatto(end));
                rho_bottom = max(rho_zLobatto(1),rho_zLobatto(end));
                dRho = rho_bottom-rho_top;
                Lz = abs(zLobatto(1)-zLobatto(end));
                
                % Now zero out the all the regions where there are density
                % inversions, project onto cheb basis, the integrate to get
                % a new density function, but keep the density at the
                % surface the same (because we use rho0 elsewhere).
                rhoz_zLobatto(rhoz_zLobatto >= 0) = max(rhoz_zLobatto(rhoz_zLobatto<0));
                rhoz_zCheb = InternalModesSpectral.fct(rhoz_zLobatto);
                rho_zCheb = (Lz/2)*InternalModesSpectral.IntegrateChebyshevVector(rhoz_zCheb);
                rho_zCheb(end) = [];
                rho_zLobatto = InternalModesSpectral.ifct(rho_zCheb);
                delta = - min(rho_zLobatto) + rho_top;
                rho_zLobatto = rho_zLobatto + delta;
                rho_zCheb(1) = rho_zCheb(1) + delta;
                
                % Re-derive all the properties
%                 new_rhoz_zLobatto = InternalModesSpectral.ifct((2/Lz)*InternalModesSpectral.DifferentiateChebyshevVector(rho_zCheb));
                
                
                
                dRho_new = abs(rho_zLobatto(end)-rho_zLobatto(1));
                dTotalVariation = (dRho_new-dRho)/dRho;
                
                if dTotalVariation < 1e-2
                    flag = 1;
                else
                    flag = 2;
                end
            end
            
            
        end
        
        function y = fInverseBisection(f, x, yMin,yMax, tol)
            %FINVERSEBISECTION(F, X)   Compute F^{-1}(X) using Bisection.
            % Taken from cumsum as part of chebfun.
            % chebfun/inv.m
            %
            % Copyright 2017 by The University of Oxford and The Chebfun Developers.
            % See http://www.chebfun.org/ for Chebfun information.
            
            a = yMin*ones(size(x));
            b = yMax*ones(size(x));
            c = (a + b)/2;
            
            while ( norm(b - a, inf) >= tol )
                vals = feval(f, c);
                % Bisection:
                I1 = ((vals-x) <= -tol);
                I2 = ((vals-x) >= tol);
                I3 = ~I1 & ~I2;
                a = I1.*c + I2.*a + I3.*c;
                b = I1.*b + I2.*c + I3.*c;
                c = (a+b)/2;
            end
            
            y = c;
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Chebyshev Methods
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [zLobatto, rho_zCheb] = ProjectOntoChebyshevPolynomialsWithTolerance(zIn, rhoFunc, tol)
            m = 3;
            m_max = 15;
            
            zMax = max(zIn);
            zMin = min(zIn);
            Lz = zMax - zMin;
                        
            n = 2^m + 1;
            cutoff = n;
            while (cutoff == n && m <= m_max)
                m = m + 1;
                n = 2^m + 1;
                
                zLobatto = (Lz/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + zMin;
                rho_zCheb = InternalModesSpectral.fct(rhoFunc(zLobatto));
                cutoff = InternalModesSpectral.standardChop(rho_zCheb, tol);
            end
            
            if cutoff < n
                n = cutoff;
                zLobatto = (Lz/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + zMin;
            else
                disp('Unable to project density function within requested tolerance! Using maximum allowed length.');
            end
            
            zLobatto(1) = zMax;
            zLobatto(end) = zMin;
            rho_zCheb = rho_zCheb(1:n);
        end
                
        % Fast Chebyshev Transform
        % By Allan P. Engsig-Karup, apek@imm.dtu.dk.
        function uh = fct(u)
            N  = length(u);
            u  = ifft(u([1:N N-1:-1:2])); % reverse ordering due to Matlab's fft
            uh = ([u(1); 2*u(2:(N-1)); u(N)]);
            
            if any(imag(uh))
                disp('Fast Chebyshev Transform returned imaginary values. Something went wrong!')
            end
            
            uh = real(uh);
        end
        
        % Inverse Fast Chebyshev Transform
        function u = ifct(uh)
            N = length(uh) - 1;
            s = N*[uh(1)*2; uh(2:N); uh(end)*2];
            u = ifft([s; flip(s(2:end-1))],'symmetric');
            u=u(1:N+1);
        end
        
        function bool = IsChebyshevGrid(z_in)
            % make sure the grid is monotonically decreasing
            if (z_in(2) - z_in(1)) > 0
                z_in = flip(z_in);
            end
            
            z_norm = InternalModesSpectral.ChebyshevPolynomialsOnGrid( z_in );
            N_points = length(z_in);
            xi=(0:N_points-1)';
            lobatto_grid = cos(xi*pi/(N_points-1));
            z_diff = z_norm-lobatto_grid;
            if max(abs(z_diff)) < 1e-6
                bool = 1;
            else
                bool = 0;
            end
        end
        
        % Given some Lobatto grid and some desired output grid, return the
        % transformation function T that goes from spectral to the output
        % grid. This basically gives you spectral interpolation.
        function [T, doesOutputGridSpanDomain] = ChebyshevTransformForGrid(lobatto_grid, output_grid)
            if(min(output_grid) < min(lobatto_grid) || max(output_grid) > max(lobatto_grid))
               error('The output grid must be bounded by the lobatto grid'); 
            end
            if (min(output_grid) == min(lobatto_grid) && max(output_grid) == max(lobatto_grid))
                doesOutputGridSpanDomain = 1;
            else
                doesOutputGridSpanDomain = 0;
            end
            
            % T_out transforms vector solutions of the eigenvalue problem
            % into gridded solution on z_out
            if doesOutputGridSpanDomain == 1 && InternalModesSpectral.IsChebyshevGrid(output_grid) == 1
                if length(output_grid) == length(lobatto_grid)
                    T = @(f_cheb) InternalModesSpectral.ifct(f_cheb);
                elseif length(output_grid) > length(lobatto_grid)
                    T = @(f_cheb) InternalModesSpectral.ifct(cat(1,f_cheb,zeros(length(output_grid)-length(lobatto_grid),1)));
                elseif length(output_grid) < length(lobatto_grid)
                    T = @(f_cheb) InternalModesSpectral.ifct(f_cheb(1:length(output_grid)));
                end
            else
                L = max(lobatto_grid)-min(lobatto_grid);
                x = (2/L)*(output_grid-min(lobatto_grid)) - 1;
                t = acos(x);
                TT = zeros(length(output_grid),length(lobatto_grid));
                for iPoly=0:(length(lobatto_grid)-1)
                    TT(:,iPoly+1) = cos(iPoly*t);
                end
                T = @(f_cheb) TT*f_cheb;
            end
        end
        
        function v_p = DifferentiateChebyshevVector(v)
            v_p = zeros(size(v));
            k = length(v)-1;
            v_p(k) = 2*k*v(k+1);
            for k=(length(v)-2):-1:1
                v_p(k) = 2*k*v(k+1) + v_p(k+2);
            end
            v_p(1) = v_p(1)/2;
        end
        
        function s = IntegrateChebyshevVectorWithLimits(v,x,a,b)
            % v are the coefficients of the chebyshev polynomials
            % x is the domain, a Gauss-Lobatto grid
            % a and b are the lower and upper limits, respectively
            v_p = InternalModesSpectral.IntegrateChebyshevVector(v);
            
            s = InternalModesSpectral.ValueOfFunctionAtPointOnGrid(b,x,v_p) - InternalModesSpectral.ValueOfFunctionAtPointOnGrid(a,x,v_p);
        end
                   
        function v_p = IntegrateChebyshevVector(v)
            % Taken from cumsum as part of chebfun.
            % chebfun/@chebtech/cumsum.m
            %
            % Copyright 2017 by The University of Oxford and The Chebfun Developers.
            % See http://www.chebfun.org/ for Chebfun information.
            
            % integration target
            n = length(v);
            v_p = zeros(n+1,1);
            
            % zero-pad coefficients
            v = reshape(v,[],1);
            v = [v; zeros(2,1)];
            
            % Compute b_(2) ... b_(n+1):
            v_p(3:n+1,:) = (v(2:n,:) - v(4:n+2,:)) ./ (2*(2:n).');
            v_p(2,:) = v(1,:) - v(3,:)/2;        % Compute b_1
            t = ones(1, n);
            t(2:2:end) = -1;
            v_p(1,:) = t*v_p(2:end,:);             % Compute b_0 (satisfies f(-1) = 0)
        end
        
        function D = ChebyshevDifferentiationMatrix(n)
            %% Chebyshev Differentiation Matrix
            % Returns the Chebyshev differentiation matrix for the first n polynomials.
            D = zeros(n,n);
            for i=1:n
                for j=1:n
                    if ( j >= i+1 && mod(i+j,2)==1 )
                        D(i,j) = 2*(j-1);
                    else
                        D(i,j) = 0.0;
                    end
                end
            end
            D(1,:)=0.5*D(1,:);
        end
        
        function D = ChebyshevInterpolationDerivative(n)
            %% Chebyshev Interpolation Derivative
            % taken from Canuto, et al. 2.4.33
            D = zeros(n,n);
            N = n-1;
            c = @(j) (j == 0 || j == N)*2 + (j>0 && j<N)*1;
            for j=0:(n-1)
                for l=0:(n-1)
                    if j ~= l
                        D(j+1,l+1) = -(c(j)/c(l))*((-1)^(j+l))/( sin( (j+l)*pi/(2*N) ) * sin( (j-l)*pi/(2*N) ) )/2;
                    elseif j == l && j == 0
                        D(j+1,l+1) = (2*N*N+1)/6;
                    elseif j == l && j == N
                        D(j+1,l+1) = -(2*N*N+1)/6;
                    else
                        D(j+1,l+1) = -cos(pi*j/N)/(2*(sin(j*pi/N))^2);
                    end
                end
            end
            D(1,:)=D(1,:);
        end
        
        function value = ValueOfFunctionAtPointOnGrid( x0, x, func_cheb )
           % We have the Chebyshev coefficents of function func_cheb, defined on grid x, return the value at x0;
           L = max(x)-min(x);
           x_norm = (2/L)*(x0-min(x)) - 1;
           t = acos(x_norm);
           
           N_polys = length(func_cheb);
%            value=sum(func_cheb.*cos(t*(0:(N_polys-1))));   
           value = func_cheb(1)*ones(size(x0));
           for i=2:N_polys
               value = value + func_cheb(i)*cos(t*(i-1));
           end
                   
        end
        
        function [varargout] = ChebyshevPolynomialsOnGrid( x, N_polys )
            %% Chebyshev Polynomials on Grid
            % Compute the the first N Chebyshev polynomials and their derivatives for
            % an arbitrary grid x.
            %
            % x_norm = ChebyshevPolynomialsOnGrid( x ) with exactly one argument, x,
            % returns the x normalized to its typical [-1,1] values.
            %
            % T = ChebyshevPolynomialsOnGrid( x, N_polys ) returns the first N_poly
            % Chebyshev polynomials for an arbitrary grid x.
            %
            % [T, T_x, T_xx,...] = ChebyshevPolynomialsOnGrid( x, N_polys ) returns the
            % first N_poly Chebyshev polynomials and their derivatives for an arbitrary
            % grid x.
            %
            % The returned matrices T, T_xx, etc are size(T) = [length(x) N_polys],
            % i.e., the polynomials are given column-wise.
            
            N_points = length(x);
            
            % These are the normalized coordinates for Chebyshev polynomials.
            L = max(x)-min(x);
            x_norm = (2/L)*(x-min(x)) - 1;
            t = acos(x_norm);
            
            % if there's only one input argument, we just return x_norm
            if nargin == 1
                varargout{1} = x_norm;
                return;
            else
                if N_polys < 4
                    disp('You must request at least four polynomials. Fixing that for you.')
                    N_polys = 4;
                end
            end
            
            N_diff = nargout-1;
            varargout = cell(1,nargout);
            
            % It's easy to create the polynomials, they're stretched cosines!
            T = zeros(N_points,N_polys);
            for iPoly=0:(N_polys-1)
                T(:,iPoly+1) = cos(iPoly*t);
            end
            
            varargout{1} = T;
            
            % Now use the recursion formula to compute derivates of polynomials.
            for n=1:N_diff
                T = varargout{n};
                T_x = zeros(size(T));
                T_x(:,2) = T(:,1);
                T_x(:,3) = 2*2*T(:,2);
                for j=4:N_polys
                    m = j-1;
                    T_x(:,j) = (m/(m-2))*T_x(:,j-2) + 2*m*T(:,j-1);
                end
                varargout{n+1} = (2/L)*T_x;
            end
            
        end
        
        function [z_lobatto_grid] = FindSmallestChebyshevGridWithNoGaps(z)
            % Want to create a chebyshev grid that never has two or more point between
            % its points. If that makes sense.
            if (z(2) - z(1)) > 0 % make z_out decreasing
                z = flip(z);
            end
            
            L = max(z)-min(z);
            np = ceil(log2(length(z)));
            n = 2^np;
            z_lobatto_grid = (L/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + min(z);
            
            while( length(unique(interp1(z_lobatto_grid,z_lobatto_grid,z,'previous'))) ~= length(z) )
                np = np + 1;
                n = 2^np;
                z_lobatto_grid = (L/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + min(z);
            end
        end
        
        function cutoff = standardChop(coeffs, tol)
            % Copyright 2017 by The University of Oxford and The Chebfun Developers. 
            % See http://www.chebfun.org/ for Chebfun information.
            %
            % This is taken, without comments and safe checks, from the
            % above developers. They get sole credit. Jared Aurentz and Nick Trefethen, July 2015.
            n = length(coeffs);
            cutoff = n;
            if ( n < 17 )
                return
            end
            
            envelope = cummax(abs(coeffs),'reverse');
            if envelope(1) == 0
                cutoff = 1;
                return
            else
                envelope = envelope/envelope(1);
            end
            
            for j = 2:n
                j2 = round(1.25*j + 5);
                if ( j2 > n )
                    % there is no plateau: exit
                    return
                end
                e1 = envelope(j);
                e2 = envelope(j2);
                r = 3*(1 - log(e1)/log(tol));
                plateau = (e1 == 0) | (e2/e1 > r);
                if ( plateau )
                    % a plateau has been found: go to Step 3
                    plateauPoint = j - 1;
                    break
                end
            end
            
            
            if ( envelope(plateauPoint) == 0 )
                cutoff = plateauPoint;
            else
                j3 = sum(envelope >= tol^(7/6));
                if ( j3 < j2 )
                    j2 = j3 + 1;
                    envelope(j2) = tol^(7/6);
                end
                cc = log10(envelope(1:j2));
                cc = cc(:);
                cc = cc + linspace(0, (-1/3)*log10(tol), j2)';
                [~, d] = min(cc);
                cutoff = max(d - 1, 1);
            end
            
        end
        
        function roots = FindRootsFromChebyshevVector(f_cheb, zLobatto)
            % Copyright (c) 2007, Stephen Morris 
            % All rights reserved.
            n = length(f_cheb);
            
            f_cheb(abs(f_cheb)<1e-15) = 1e-15;
            
            A=zeros(n-1);   % "n-1" because Boyd's notation includes the zero-indexed
            A(1,2)=1;       % elements whereas Matlab's of course does not allow this.
            % In other words, to replicate Boyd's N=5 matrix we need to
            % set n=6.
            for j=2:n-2
                for k=1:n-1
                    if j==k+1 || j==k-1
                        A(j,k)=0.5;
                    end
                end
            end
            for k=1:n-1
                A(n-1,k)=-f_cheb(k)/(2*f_cheb(n));  % c(1) in our notation is c(0) in Boyd's
            end
            A(n-1,n-2)=A(n-1,n-2)+0.5;
            % Now we have the companion matrix, we can find its eigenvalues using the
            % MATLAB built-in function.
            eigvals=eig(A);
            
            % We're only interested in the real elements of the matrix:
            realvals=(arrayfun(@(x) ~any(imag(x)),eigvals)).*eigvals;
            
            % Of course these are the roots scaled to the canonical interval [-1,1]. We
            % need to map them back onto the interval [a,b]; we widen the interval just
            % a fraction to make sure that we don't miss any that are right on the
            % edge.
            a = min(zLobatto);
            b = max(zLobatto);
            rangevals=nonzeros((arrayfun(@(x) abs(x)<=1.001, realvals)).*realvals);
            roots=sort((rangevals.*0.5*(b-a)) + (0.5*(b+a)));
        end
        
    end
end


