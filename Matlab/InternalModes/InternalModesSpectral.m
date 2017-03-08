classdef InternalModesSpectral < InternalModesBase
    properties (Access = public)
        rho
        N2
    end
    
    properties (Dependent)
        rho_z
        rho_zz
    end
    
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
        
        function f = SetNoiseFloorToZero(~, f)
           f(abs(f)/max(abs(f)) < 1e-15) = 0;
        end
        
        % Called by InitializeWithGrid and InitializeWithFunction after
        % they've completed their basic setup.
        function self = InitializeDerivedQuanitites(self)
            self.rho_zCheb = fct(self.rho_zLobatto);
            self.rho_zCheb = self.SetNoiseFloorToZero(self.rho_zCheb);
            self.Diff1_zCheb = (2/self.Lz)*ChebyshevDifferentiationMatrix( length(self.zLobatto) );
            self.N2_zCheb= -(self.g/self.rho0)*self.Diff1_zCheb*self.rho_zCheb;
            
            N2_zLobatto = ifct(self.N2_zCheb);
            if any(N2_zLobatto < 0)
               error('The bouyancy frequency goes negative! This is likely happening because of spline interpolation of density. Try using a finer grid.'); 
            end
            
            self.T_zCheb_zOut = ChebyshevTransformForGrid(self.zLobatto, self.z);
            self.rho = self.T_zCheb_zOut(self.rho_zCheb);
            self.N2 = self.T_zCheb_zOut(self.N2_zCheb);
        end
        
        % Superclass calls this method upon initialization when it
        % determines that the input is given in gridded form. The goal is
        % to initialize zLobatto and rho_zLobatto.
        function self = InitializeWithGrid(self, rho, zIn)
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
            if (length(self.zLobatto) < 1024)
                n = 1024;
                self.zLobatto = (self.Lz/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + zMin;
            end
            
            self.rho_zLobatto = rho(self.zLobatto);
            
            self.InitializeDerivedQuanitites();
        end
        
        function self = SetupEigenvalueProblem(self)
            % The user requested that the eigenvalue problem be solved on a
            % grid of particular length
            if self.nEVP > 0
                if self.nModes > self.nEVP
                    self.nEVP = self.nModes;
                end
            else
                self.nEVP = 512;
            end
            
            n = self.nEVP;
            self.xLobatto = (self.Lz/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + min(self.zLobatto);
            T_zCheb_xLobatto = ChebyshevTransformForGrid(self.zLobatto, self.xLobatto);
            
            self.N2_xLobatto = T_zCheb_xLobatto(self.N2_zCheb);
            self.Diff1_xCheb = (2/self.Lz)*ChebyshevDifferentiationMatrix( length(self.xLobatto) );
            
            [self.T_xLobatto,self.Tx_xLobatto,self.Txx_xLobatto] = ChebyshevPolynomialsOnGrid( self.xLobatto, length(self.xLobatto) );
            [self.T_xCheb_zOut, self.doesOutputGridSpanDomain] = ChebyshevTransformForGrid(self.xLobatto, self.z);
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
            
            hFromLambda = @(lambda) 1.0 ./ lambda;
            GFromGCheb = @(G_cheb,h) self.T_xCheb_zOut(G_cheb);
            FFromGCheb = @(G_cheb,h) h * self.T_xCheb_zOut(self.Diff1_xCheb*G_cheb);
            [F,G,h] = ModesFromGEPSpecial(self,A,B,hFromLambda,GFromGCheb,FFromGCheb);
        end
        
        function [F,G,h] = ModesAtFrequency(self, omega )
            error('This function is not yet implemented!');
        end
        
        % Take matrices A and B from the generalized eigenvalue problem
        % (GEP) and returns F,G,h. The h_func parameter is a function that
        % returns the eigendepth, h, given eigenvalue lambda from the GEP.
        function [F,G,h] = ModesFromGEP(self,A,B,hFromLambda,GFromGCheb,FFromGCheb)
            [V,D] = eig( A, B );
            
            [lambda, permutation] = sort(abs(diag(D)),1,'ascend');
            G_cheb=V(:,permutation);
            h = hFromLambda(lambda.');
            
            n = size(G_cheb,1);
            F = zeros(length(self.z),n);
            G = zeros(length(self.z),n);
            
            if self.doesOutputGridSpanDomain == 1
                for j=1:n
                    G(:,j) = GFromGCheb(G_cheb(:,j),h(j));
                    F(:,j) = FFromGCheb(G_cheb(:,j),h(j));
                end
                [F,G] = self.NormalizeModes(F,G,self.N2,self.z);
            else
                error('This normalization condition is not yet implemented!')
            end
        end
        
        % In theory we can integrate spectrally, but this seems to do worse
        % than trapz
        function [F,G,h] = ModesFromGEPSpecial(self,A,B,hFromLambda,GFromGCheb,FFromGCheb)
            [V,D] = eig( A, B );
            
            [lambda, permutation] = sort(abs(diag(D)),1,'ascend');
            G_cheb=V(:,permutation);
            h = hFromLambda(lambda.');
            
            n = size(G_cheb,1);
            F = zeros(length(self.z),n);
            G = zeros(length(self.z),n);
            
            np = (0:(n-1))';
            Int1 = -(1+(-1).^np)./(np.*np-1);
            Int1(2) = 0;
            Int1 = self.Lz/2*Int1;
            
            [~,maxIndexZ] = max(self.z);
            if self.doesOutputGridSpanDomain == 1
                for j=1:n
                    J = (1/self.g) * (self.N2_xLobatto - self.f0*self.f0) .* ifct(G_cheb(:,j)) .^ 2;
                    A = sqrt(abs(sum(Int1.*fct(J))));
                    
                    G(:,j) = GFromGCheb(G_cheb(:,j),h(j))/A;
                    F(:,j) = FFromGCheb(G_cheb(:,j),h(j))/A;
                    
                    if F(maxIndexZ,j)< 0
                        F(:,j) = -F(:,j);
                        G(:,j) = -G(:,j);
                    end
                end
%                 [F,G] = self.NormalizeModes(F,G,self.z);
            else
                error('This normalization condition is not yet implemented!')
            end
        end
        
         
    end
    
end


