classdef InternalModesSpectral < InternalModesBase
    % Need to establish notation:
    %
    % Dimensions ? like z ? have grids, which may differ. zIn is the
    % dimension z, on some grid given by the user. zLobatto is the
    % diension z, on a Lobatto grid. zCheb would be used for the Chebyshev
    % polynomial representation of a function defined on z. Note that the
    % base class property 'z' is really zOut.
    %
    % Functions ? if they're analytical, they're defined on a dimension,
    % not a grid. If they're gridded, then the grid must also be specified.
    % Thus, the function rho(z) when analytically defined will be given as
    % either just rho (when obvious) or rho_z when necessary. If it's
    % gridded, then expect rho_zIn or rho_zLobatto.
    %
    % Transformations ? a transformation often a matrix (or generally just
    % a function) that transformation a function from one basis to another.
    % A common scenario here is that we need to go from zLobatto to zCheb.
    % In our cases we're usually going from one gridded dimension to
    % another. Inverse Chebyshev transformations are denote by T, or more
    % specifically T_zCheb_zLobatto. Because T denotes Chebyshev, we will
    % simply reduce this to T_zLobatto
    %
    % The underscore in the case of rho_z is certainly recognized as a
    % derivative (because of LaTex notation), so it is awkward with the
    % notion used here, where the underscore is used to separate the
    % function/transformation name from the grid. Thus, we leave rho_z and
    % rho_zz as the only two exceptions to the aforementioned notation
    % because they are the only public facing interface.
    properties (Access = public)
        rho
        N2
    end
    
    properties (Dependent)
        rho_z
        rho_zz
    end
    

    properties %(Access = private)
        zLobatto            % zIn coordinate, on Chebyshev extrema/Lobatto grid
        rho_zLobatto        % rho on the above zLobatto
        rho_zCheb           % rho in spectral space
        N2_zCheb            % N2 in spectral space
        Diff1_zCheb         % single derivative in spectral space
        T_zCheb_zOut        % *function* handle that transforms from spectral to z_out
        
        % These are the grids used for solving the eigenvalue problem.
        % Really x=z, although they may have different numbers of points.
        nEVP = 0           % number of points in the eigenvalue problem
        xLobatto           % Lobatto grid on z with nEVP points. May be the same as zLobatto, may not be.
        N2_xLobatto        % N2 on the z2Lobatto grid
        Diff1_xCheb        % single derivative in spectral space
        T_xLobatto, Tx_xLobatto, Txx_xLobatto        % Chebyshev polys (and derivs) on the zLobatto
        T_xCheb_zOut
        Int_xCheb           % Vector that multiplies Cheb coeffs, then sum for integral
        doesOutputGridSpanDomain = 0    % if not, the output grid can't be used for normalization
        
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
            value = self.T_zCheb_zOut(self.Diff1_zCheb * self.rho_zCheb);
        end
        
        function value = get.rho_zz(self)
            value = self.T_zCheb_zOut(self.Diff1_zCheb * self.Diff1_zCheb * self.rho_zCheb);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computation of the modes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [F,G,h] = ModesAtWavenumber(self, k )
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
        
        function f = SetNoiseFloorToZero(~, f)
            f(abs(f)/max(abs(f)) < 1e-15) = 0;
        end
        
        % Called by InitializeWithGrid and InitializeWithFunction after
        % they've completed their basic setup.
        function self = InitializeDerivedQuanitites(self)
            self.rho_zCheb = InternalModesSpectral.fct(self.rho_zLobatto);
            self.rho_zCheb = self.SetNoiseFloorToZero(self.rho_zCheb);
            self.Diff1_zCheb = (2/self.Lz)*InternalModesSpectral.ChebyshevDifferentiationMatrix( length(self.zLobatto) );
            self.N2_zCheb= -(self.g/self.rho0)*self.Diff1_zCheb*self.rho_zCheb;
            
            N2_zLobatto = InternalModesSpectral.ifct(self.N2_zCheb);
            if any(N2_zLobatto < 0)
                fprintf('The bouyancy frequency goes negative! This is likely happening because of spline interpolation of density. We will proceed by setting N2=0 at those points.\n');
                N2_zLobatto(N2_zLobatto < 0) = 0;
                self.N2_zCheb = InternalModesSpectral.fct(N2_zLobatto);
            end
            
            self.T_zCheb_zOut = InternalModesSpectral.ChebyshevTransformForGrid(self.zLobatto, self.z);
            self.rho = self.T_zCheb_zOut(self.rho_zCheb);
            self.N2 = self.T_zCheb_zOut(self.N2_zCheb);
        end
        
        % Superclass calls this method upon initialization when it
        % determines that the input is given in gridded form. The goal is
        % to initialize zLobatto and rho_zLobatto.
        function self = InitializeWithGrid(self, rho, zIn)
            self.validateInitialModeAndEVPSettings();
            
            if (zIn(2) - zIn(1)) > 0
                zIn = flip(zIn);
                rho = flip(rho);
            end
            
            if IsChebyshevGrid(zIn) == 1
                self.zLobatto = zIn;
                self.rho_zLobatto = rho;
            else
                self.zLobatto = FindSmallestChebyshevGridWithNoGaps( zIn ); % z, on a chebyshev grid
                self.rho_zLobatto = interp1(zIn, rho, self.zLobatto, 'spline'); % rho, interpolated to that grid
            end
            
            self.InitializeDerivedQuanitites();
        end
        
        % Superclass calls this method upon initialization when it
        % determines that the input is given in functional form. The goal
        % is to initialize zLobatto and rho_zLobatto.
        function self = InitializeWithFunction(self, rho, zMin, zMax, zOut)
            
            self.validateInitialModeAndEVPSettings();
            
            % We have an analytical function describing density, so lets
            % place points on a Chebyshev grid such that z_out is fully
            % resolved.
            if (zOut(2) - zOut(1)) > 0 % make z_out decreasing
                zOut = flip(zOut);
            end
            
            zIn = zOut;
            % Note that z_out might not span the whole domain, so we may
            % need to add z_min and z_max to the end points.
            if zMax > zIn(1)
                zIn = cat(1,zMax,zIn);
            end
            if zMin < zIn(end)
                zIn = cat(1,zIn,zMin);
            end
            
            if IsChebyshevGrid(zIn) == 1
                self.zLobatto = zIn;
            else
                self.zLobatto = FindSmallestChebyshevGridWithNoGaps( zIn ); % z, on a chebyshev grid
            end
            
            % There's just no reason not to go at least this big since
            % this is all fast transforms.
            if (length(self.zLobatto) < 1025)
                n = 1025; % 2^n + 1 for a fast Chebyshev transform
                self.zLobatto = (self.Lz/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + zMin;
            end
            
            self.rho_zLobatto = rho(self.zLobatto);
            
            self.InitializeDerivedQuanitites();
        end
        
        function self = SetupEigenvalueProblem(self)
            n = self.nEVP;
            self.xLobatto = (self.Lz/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + min(self.zLobatto);
            T_zCheb_xLobatto = InternalModesSpectral.ChebyshevTransformForGrid(self.zLobatto, self.xLobatto);
            
            self.N2_xLobatto = T_zCheb_xLobatto(self.N2_zCheb);
            self.Diff1_xCheb = (2/self.Lz)*InternalModesSpectral.ChebyshevDifferentiationMatrix( length(self.xLobatto) );
            
            [self.T_xLobatto,self.Tx_xLobatto,self.Txx_xLobatto] = InternalModesSpectral.ChebyshevPolynomialsOnGrid( self.xLobatto, length(self.xLobatto) );
            [self.T_xCheb_zOut, self.doesOutputGridSpanDomain] = InternalModesSpectral.ChebyshevTransformForGrid(self.xLobatto, self.z);
            
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
            
            if isempty(self.nModes) || self.nModes < 1
                self.nModes = floor(self.nEVP/2);
            end
        end
        
        % This function is an intermediary used by ModesAtFrequency and
        % ModesAtWavenumber to establish the various norm functions.
        function [F,G,h] = ModesFromGEPSpectral(self,A,B)
            hFromLambda = @(lambda) 1.0 ./ lambda;
            GOutFromGCheb = @(G_cheb,h) self.T_xCheb_zOut(G_cheb);
            FOutFromGCheb = @(G_cheb,h) h * self.T_xCheb_zOut(self.Diff1_xCheb*G_cheb);
            GFromGCheb = @(G_cheb,h) InternalModesSpectral.ifct(G_cheb);
            FFromGCheb = @(G_cheb,h) h * InternalModesSpectral.ifct( self.Diff1_xCheb*G_cheb );
            GNorm = @(Gj) abs(sum(self.Int_xCheb .*InternalModesSpectral.fct((1/self.g) * (self.N2_xLobatto - self.f0*self.f0) .* Gj .^ 2)));
            FNorm = @(Fj) abs(sum(self.Int_xCheb .*InternalModesSpectral.fct((1/self.Lz) * Fj.^ 2)));
            [F,G,h] = ModesFromGEP(self,A,B,hFromLambda,GFromGCheb,FFromGCheb,GNorm,FNorm,GOutFromGCheb,FOutFromGCheb);
        end
        
        % Take matrices A and B from the generalized eigenvalue problem
        % (GEP) and returns F,G,h. The h_func parameter is a function that
        % returns the eigendepth, h, given eigenvalue lambda from the GEP.
        function [F,G,h] = ModesFromGEP(self,A,B,hFromLambda,GFromGCheb, FFromGCheb, GNorm,FNorm, GOutFromGCheb,FOutFromGCheb)
            [V,D] = eig( A, B );
            
            [lambda, permutation] = sort(abs(diag(D)),1,'ascend');
            G_cheb=V(:,permutation);
            h = hFromLambda(lambda.');
            
            F = zeros(length(self.z),self.nModes);
            G = zeros(length(self.z),self.nModes);
            h = h(1:self.nModes);
            
            % This still need to be optimized to *not* do the transforms
            % twice, when the EVP grid is the same as the output grid.
            [~,maxIndexZ] = max(self.zLobatto);
            for j=1:self.nModes
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
        
        % Given some Lobatto grid and some desired output grid, return the
        % transformation function T that goes from spectral to the output
        % grid. This basically gives you spectral interpolation.
        function [T, doesOutputGridSpanDomain] = ChebyshevTransformForGrid(lobatto_grid, output_grid)
            if (min(output_grid) == min(lobatto_grid) && max(output_grid) == max(lobatto_grid))
                doesOutputGridSpanDomain = 1;
            else
                doesOutputGridSpanDomain = 0;
            end
            
            % T_out transforms vector solutions of the eigenvalue problem
            % into gridded solution on z_out
            if doesOutputGridSpanDomain == 1 && IsChebyshevGrid(output_grid) == 1
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
    end
end


