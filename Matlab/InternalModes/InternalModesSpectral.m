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
        rho % Density on the z grid.
        N2 % Buoyancy frequency on the z grid, $N^2 = -\frac{g}{\rho(0)} \frac{\partial \rho}{\partial z}$.
    end
    
    properties (Dependent)
        rho_z % First derivative of density on the z grid.
        rho_zz % Second derivative of density on the z grid.
    end
    

    properties %(Access = private)
        rho_function        % function handle to return rho at z (uses interpolation for grids)
        
        % These properties are initialized with InitializeZLobattoProperties()
        % Most subclasses would *not* override these initializations.
        zLobatto            % zIn coordinate, on Chebyshev extrema/Lobatto grid
        rho_zLobatto        % rho on the above zLobatto
        rho_zCheb           % rho in spectral space
        N2_zCheb            % N2 in spectral space
        N2_zLobatto         % N2 on the zLobatto grid
        Diff1_zCheb         % *function handle*, that takes single derivative in spectral space
        T_zCheb_zOut        % *function handle*, that transforms from spectral to z_out
        
        % These properties are initialized with InitializeXLobattoProperties()
        % Most subclasses *will* override these initializations. The 'x' refers to the stretched coordinate being used.
        % This class uses x=z (depth), although they may have different numbers of points.
        nEVP = 0           % number of points in the eigenvalue problem
        nGrid = 0          % number of points used to compute the derivatives of density.
        xLobatto           % stretched coordinate Lobatto grid nEVP points. z for this class, density or wkb for others.
        z_xLobatto         % The value of z, at the xLobatto points
        xOut               % desired locations of the output in x-coordinate (deduced from z_out)
        N2_xLobatto        % N2 on the z2Lobatto grid
        Diff1_xCheb        % single derivative in spectral space, *function handle*
        T_xLobatto, Tx_xLobatto, Txx_xLobatto        % Chebyshev polys (and derivs) on the zLobatto
        T_xCheb_zOut
        Int_xCheb           % Vector that multiplies Cheb coeffs, then sum for integral
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
            
            fprintf('N2 was computed with %d points. The eigenvalue problem will be computed with %d points.\n',length(self.zLobatto), length(self.xLobatto));
        end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computed (dependent) properties
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function value = get.rho_z(self)
            value = self.T_zCheb_zOut(self.Diff1_zCheb(self.rho_zCheb));
        end
        
        function value = get.rho_zz(self)
            value = self.T_zCheb_zOut(self.Diff1_zCheb(self.Diff1_zCheb(self.rho_zCheb)));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computation of the modes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [F,G,h] = ModesAtWavenumber(self, k )
            self.gridFrequency = 0;
            
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
            
            % Lower boundary is rigid, G=0
            A(n,:) = T(n,:);
            B(n,:) = 0;
            
            % G=0 or G_z = \frac{1}{h_j} G at the surface, depending on the BC
            if strcmp(self.upperBoundary, 'free_surface')
                % G_z = \frac{1}{h_j} G at the surface
                A(1,:) = Tz(1,:);
                B(1,:) = T(1,:);
            elseif strcmp(self.upperBoundary, 'rigid_lid')
                A(1,:) = T(1,:);
                B(1,:) = 0;
            end
            
            [F,G,h] = self.ModesFromGEPSpectral(A,B);
        end
        
        function [F,G,h] = ModesAtFrequency(self, omega )
            self.gridFrequency = omega;
            
            T = self.T_xLobatto;
            Tz = self.Tx_xLobatto;
            Tzz = self.Txx_xLobatto;
            n = self.nEVP;
            
            A = Tzz;
            B = diag((omega*omega-self.N2_xLobatto)/self.g)*T;
            
            % Lower boundary is rigid, G=0
            A(n,:) = T(n,:);
            B(n,:) = 0;
            
            % G=0 or G_z = \frac{1}{h_j} G at the surface, depending on the BC
            if strcmp(self.upperBoundary, 'free_surface')
                % G_z = \frac{1}{h_j} G at the surface
                A(1,:) = Tz(1,:);
                B(1,:) = T(1,:);
            elseif strcmp(self.upperBoundary, 'rigid_lid')
                A(1,:) = T(1,:);
                B(1,:) = 0;
            end
            
            [F,G,h] = self.ModesFromGEPSpectral(A,B);
        end
        
 
    end
    
    methods (Access = protected)       
        function self = InitializeWithGrid(self, rho, zIn)
            % Superclass calls this method upon initialization when it
            % determines that the input is given in gridded form. The goal is
            % to initialize zLobatto and rho_zLobatto.
            self.validateInitialModeAndEVPSettings();
            
            % Our Chebyshev grids are monotonically decreasing.
            if (zIn(2) - zIn(1)) > 0
                zIn = flip(zIn);
                rho = flip(rho);
            end
                      
            self.rho_function = @(z) interp1(zIn, rho, z, 'spline');
            if InternalModesSpectral.IsChebyshevGrid(zIn) == 1
                self.zLobatto = zIn;
                self.rho_zLobatto = rho;
            else
                self.zLobatto = InternalModesSpectral.FindSmallestChebyshevGridWithNoGaps( zIn ); % z, on a chebyshev grid
                self.rho_zLobatto = self.rho_function(self.zLobatto); % rho, interpolated to that grid
            end
              
            self.InitializeZLobattoProperties();
        end
        
        function self = InitializeWithFunction(self, rho, zMin, zMax, zOut)
            % Superclass calls this method upon initialization when it
            % determines that the input is given in functional form. The goal
            % is to initialize zLobatto and rho_zLobatto.
            self.validateInitialModeAndEVPSettings();
                        
            if self.nGrid < 65 || isempty(self.nGrid) 
                self.nGrid = 2^14 + 1; % 2^n + 1 for a fast Chebyshev transform
            end
            
            n = self.nGrid;
            self.zLobatto = (self.Lz/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + zMin;  
            self.rho_zLobatto = rho(self.zLobatto);
            self.rho_function = rho;
            
            self.InitializeZLobattoProperties();
        end      
        
        % Called by InitializeWithGrid and InitializeWithFunction after
        % they've completed their basic setup.
        function self = InitializeZLobattoProperties(self)
            self.rho_zCheb = InternalModesSpectral.fct(self.rho_zLobatto);
            self.rho_zCheb = self.SetNoiseFloorToZero(self.rho_zCheb);
            self.Diff1_zCheb = @(v) (2/self.Lz)*InternalModesSpectral.DifferentiateChebyshevVector( v );
            self.N2_zCheb= -(self.g/self.rho0)*self.Diff1_zCheb(self.rho_zCheb);
                        
            self.N2_zLobatto = InternalModesSpectral.ifct(self.N2_zCheb);
            if any(self.N2_zLobatto < 0)
                fprintf('The bouyancy frequency goes negative! This is likely happening because of spline interpolation of density. We will proceed by setting N2=0 at those points.\n');
                self.N2_zLobatto(self.N2_zLobatto <= 0) = min(self.N2_zLobatto(self.N2_zLobatto>0));
                self.N2_zCheb = InternalModesSpectral.fct(self.N2_zLobatto);
            end
            
            self.T_zCheb_zOut = InternalModesSpectral.ChebyshevTransformForGrid(self.zLobatto, self.z);
            self.rho = self.T_zCheb_zOut(self.rho_zCheb);
            self.N2 = self.T_zCheb_zOut(self.N2_zCheb);
        end
        
        function self = SetupEigenvalueProblem(self)
            n = self.nEVP;
            self.xLobatto = (self.Lz/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + min(self.zLobatto);
            T_zCheb_xLobatto = InternalModesSpectral.ChebyshevTransformForGrid(self.zLobatto, self.xLobatto);
            
            self.N2_xLobatto = T_zCheb_xLobatto(self.N2_zCheb);
            self.Diff1_xCheb = @(v) (2/self.Lz)*InternalModesSpectral.DifferentiateChebyshevVector( v );
            
            [self.T_xLobatto,self.Tx_xLobatto,self.Txx_xLobatto] = InternalModesSpectral.ChebyshevPolynomialsOnGrid( self.xLobatto, length(self.xLobatto) );
            self.T_xCheb_zOut = InternalModesSpectral.ChebyshevTransformForGrid(self.xLobatto, self.z);
            
            % We use that \int_{-1}^1 T_n(x) dx = \frac{(-1)^n + 1}{1-n^2}
            % for all n, except n=1, where the integral is zero.
            np = (0:(n-1))';
            self.Int_xCheb = -(1+(-1).^np)./(np.*np-1);
            self.Int_xCheb(2) = 0;
            self.Int_xCheb = self.Lz/2*self.Int_xCheb;
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
        
        function f = SetNoiseFloorToZero(~, f)
            f(abs(f)/max(abs(f)) < 1e-15) = 0;
        end
        
        function [zBoundariesAndTPs, thesign, boundaryIndices] = FindTurningPointBoundariesAtFrequency(self, omega)
            % This function returns not just the turning points, but also
            % the top and bottom boundary locations in z. The boundary
            % indices are the index to the point just *above* the turning
            % point.
            N2Omega2_zLobatto = self.N2_zLobatto - omega*omega;
            a = N2Omega2_zLobatto; a(a>=0) = 1; a(a<0) = 0;
            turningIndices = find(diff(a)~=0);
            nTP = length(turningIndices);
            zTP = zeros(nTP,1);
            for i=1:nTP
                fun = @(z) interp1(self.zLobatto,N2Omega2_zLobatto,z,'spline');
                z0 = [self.zLobatto(turningIndices(i)+1) self.zLobatto(turningIndices(i))];
                zTP(i) = fzero(fun,z0);
            end
            boundaryIndices = [1; turningIndices; length(self.zLobatto)];
            zBoundariesAndTPs = [self.zLobatto(1); zTP; self.zLobatto(end)];
            
            % what's the sign on the EVP in these regions? In this case,
            % positive indicates oscillatory, negative exponential decay
            midZ = zBoundariesAndTPs(1:end-1) + diff(zBoundariesAndTPs)/2;
            thesign = sign( interp1(self.zLobatto,N2Omega2_zLobatto,midZ,'spline') );
            % remove any boundaries of zero length
            for index = reshape(find( thesign == 0),1,[])
                thesign(index) = [];
                zBoundariesAndTPs(index) = [];
                boundaryIndices(index) = [];
            end
        end
        
        % This function is an intermediary used by ModesAtFrequency and
        % ModesAtWavenumber to establish the various norm functions.
        function [F,G,h] = ModesFromGEPSpectral(self,A,B)
            hFromLambda = @(lambda) 1.0 ./ lambda;
            GOutFromGCheb = @(G_cheb,h) self.T_xCheb_zOut(G_cheb);
            FOutFromGCheb = @(G_cheb,h) h * self.T_xCheb_zOut(self.Diff1_xCheb(G_cheb));
            GFromGCheb = @(G_cheb,h) InternalModesSpectral.ifct(G_cheb);
            FFromGCheb = @(G_cheb,h) h * InternalModesSpectral.ifct( self.Diff1_xCheb(G_cheb) );
            GNorm = @(Gj) abs(sum(self.Int_xCheb .*InternalModesSpectral.fct((1/self.g) * (self.N2_xLobatto - self.f0*self.f0) .* Gj .^ 2)));
            FNorm = @(Fj) abs(sum(self.Int_xCheb .*InternalModesSpectral.fct((1/self.Lz) * Fj.^ 2)));
            [F,G,h] = ModesFromGEP(self,A,B,hFromLambda,GFromGCheb,FFromGCheb,GNorm,FNorm,GOutFromGCheb,FOutFromGCheb);
        end
        
        % Take matrices A and B from the generalized eigenvalue problem
        % (GEP) and returns F,G,h. The h_func parameter is a function that
        % returns the eigendepth, h, given eigenvalue lambda from the GEP.
        function [F,G,h] = ModesFromGEP(self,A,B,hFromLambda,GFromGCheb, FFromGCheb, GNorm,FNorm, GOutFromGCheb,FOutFromGCheb)
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
            
            % This still need to be optimized to *not* do the transforms
            % twice, when the EVP grid is the same as the output grid.
            [maxIndexZ] = find(self.N2_xLobatto-self.gridFrequency*self.gridFrequency>0,1,'first');
            if maxIndexZ > 1 % grab a point just above the turning point, which should have the right sign.
                maxIndexZ = maxIndexZ-1;
            end
            for j=1:maxModes
                Fj = FFromGCheb(G_cheb(:,j),h(j));
                Gj = GFromGCheb(G_cheb(:,j),h(j));
                if strcmp(self.normalization, 'max_u')
                    A = max( abs( Fj ));
                elseif strcmp(self.normalization, 'max_w')
                    A = max( abs( Gj ) );
                elseif strcmp(self.normalization, 'const_G_norm')
                    A = sqrt(GNorm( Gj ));
                elseif strcmp(self.normalization, 'const_F_norm')
                    A = sqrt(FNorm( Fj ));
                end
                if Fj(maxIndexZ) < 0
                    A = -A;
                end
                
                G(:,j) = GOutFromGCheb(G_cheb(:,j),h(j))/A;
                F(:,j) = FOutFromGCheb(G_cheb(:,j),h(j))/A;
            end
        end
    end
    
    methods (Static)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Chebyshev Methods
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
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
        
        function value = ValueOfFunctionAtPointOnGrid( x0, x, func_cheb )
           % We have the Chebyshev coefficents of function func_cheb, defined on grid x, return the value at x0;
           L = max(x)-min(x);
           x_norm = (2/L)*(x0-min(x)) - 1;
           t = acos(x_norm);
           
           N_polys = length(func_cheb);
           value=sum(func_cheb.*cos(t*(0:(N_polys-1))));           
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
    end
end


