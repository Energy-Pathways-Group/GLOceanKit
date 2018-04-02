classdef InternalModesDensitySpectral < InternalModesSpectral
    % InternalModesDensitySpectral This class solves the vertical
    % eigenvalue problem on a stretched density coordinate grid using
    % Chebyshev polynomials.
    %
    % See InternalModesBase for basic usage information.
    %
    % This class uses the coordinate s=-g*rho/rho0 to solve the EVP.
    % 
    % Internally, sLobatto is the stretched density coordinate on a Chebyshev
    % extrema/Lobatto grid. This is the grid upon which the eigenvalue problem
    % is solved, and therefore the class uses the superclass properties denoted
    % with 'x' when setting up the eigenvalue problem.
    %
    %   See also INTERNALMODES, INTERNALMODESBASE, INTERNALMODESSPECTRAL,
    %   INTERNALMODESWKBSPECTRAL, and INTERNALMODESFINITEDIFFERENCE.
    %
    %   Jeffrey J. Early
    %   jeffrey@jeffreyearly.com
    %
    %   March 14th, 2017        Version 1.0
    
    properties %(Access = private)            
        N2z_xLobatto    	% (d/dz)N2 on the z_sLobatto grid   
    end

    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = InternalModesDensitySpectral(rho, z_in, z_out, latitude, varargin)
            self@InternalModesSpectral(rho,z_in,z_out,latitude, varargin{:});
        end
                    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computation of the modes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [F,G,h,omega] = ModesAtWavenumber(self, k )
            self.gridFrequency = 0;
            
            T = self.T_xLobatto;
            Tz = self.Tx_xLobatto;
            Tzz = self.Txx_xLobatto;
            n = self.nEVP;
                         
            A = diag(self.N2_xLobatto .* self.N2_xLobatto)*Tzz + diag(self.N2z_xLobatto)*Tz - k*k*T;
            B = diag( (self.f0*self.f0 - self.N2_xLobatto)/self.g )*T;
            
            iSurface = 1;
            iBottom = n;
            
            % Lower boundary is rigid, G=0
            A(iBottom,:) = T(iBottom,:);
            B(iBottom,:) = 0;
            
            % G=0 or N^2 G_s = \frac{1}{h_j} G at the surface, depending on the BC
            if self.upperBoundary == UpperBoundary.freeSurface
                % G_z = \frac{1}{h_j} G at the surface
                A(iSurface,:) = self.N2_xLobatto(iSurface) * Tz(iSurface,:);
                B(iSurface,:) = T(iSurface,:);
            elseif self.upperBoundary == UpperBoundary.rigidLid
                A(iSurface,:) = T(iSurface,:);
                B(iSurface,:) = 0;
            end
            
            [F,G,h] = self.ModesFromGEPDensitySpectral(A,B);
            omega = self.omegaFromK(h,k);
        end
        
        function [F,G,h,k] = ModesAtFrequency(self, omega )
            self.gridFrequency = omega;
            
            T = self.T_xLobatto;
            Tz = self.Tx_xLobatto;
            Tzz = self.Txx_xLobatto;
            n = self.nEVP;
                         
            A = diag(self.N2_xLobatto .* self.N2_xLobatto)*Tzz + diag(self.N2z_xLobatto)*Tz;
            B = diag( (omega*omega - self.N2_xLobatto)/self.g )*T;
            
            iSurface = 1;
            iBottom = n;
            
            % Lower boundary is rigid, G=0
            A(iBottom,:) = T(iBottom,:);
            B(iBottom,:) = 0;
            
            % G=0 or N^2 G_s = \frac{1}{h_j} G at the surface, depending on the BC
            if self.upperBoundary == UpperBoundary.freeSurface
                % G_z = \frac{1}{h_j} G at the surface
                A(iSurface,:) = self.N2_xLobatto(iSurface) * Tz(iSurface,:);
                B(iSurface,:) = T(iSurface,:);
            elseif self.upperBoundary == UpperBoundary.rigidLid
                A(iSurface,:) = T(iSurface,:);
                B(iSurface,:) = 0;
            end
            
            [F,G,h] = self.ModesFromGEPDensitySpectral(A,B);
            k = self.kFromOmega(h,omega);
        end

    end
    
    methods (Access = protected)
        % Superclass calls this method upon initialization when it
        % determines that the input is given in gridded form.
        %
        % The superclass will initialize zLobatto and rho_lobatto;
        % this class must initialize the sLobatto, z_sLobatto and
        % sOut.
        function self = InitializeWithGrid(self, rho, z_in)
            % Superclass initializes zLobatto and rho_lobatto
            InitializeWithGrid@InternalModesSpectral(self, rho, z_in);
                        
            n = self.nEVP;
            s = -self.g*rho/self.rho0 + self.g;
            Ls = max(s)-min(s);
            self.xLobatto = (Ls/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + min(s);
            self.z_xLobatto = interp1(s, z_in, self.xLobatto, 'spline'); % z, evaluated on that s grid
            self.xOut = interp1(z_in, s, self.z, 'spline'); % z, evaluated on that s grid
        end
        
        % Superclass calls this method upon initialization when it
        % determines that the input is given in functional form.
        %
        % The superclass will initialize zLobatto and rho_lobatto;
        % this class must initialize the sLobatto, z_sLobatto and
        % sOut.
        function self = InitializeWithFunction(self, rho, zMin, zMax, zOut)
            % Superclass initializes zLobatto and rho_lobatto
            InitializeWithFunction@InternalModesSpectral(self, rho, zMin, zMax, zOut);
            
            % Create a stretched grid that includes all the points of z_out
            s = @(z) -self.g*rho(z)/self.rho0 + self.g;   
            
            n = self.nEVP;
            sMin = min([s(zMax), s(zMin)]);
            sMax = max([s(zMax), s(zMin)]);
            Ls = sMax-sMin;
            self.xLobatto = (Ls/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + sMin;
            
            % Now create a transformation for functions defined on
            % z_lobatto to (spectrally) take them into s_lobatto.
            % We use the fact that we have a function handle to iteratively
            % improve this projection.
            self.z_xLobatto = interp1(s(self.zLobatto), self.zLobatto, self.xLobatto, 'linear', self.zLobatto(1));
            nloops = 0;
            while ( nloops < 10 && max( abs(s(self.z_xLobatto) - self.xLobatto))/max(abs(self.xLobatto)) > 1e-15)
                self.z_xLobatto = interp1(s(self.z_xLobatto), self.z_xLobatto, self.xLobatto, 'linear', self.zLobatto(1));
                nloops = nloops + 1;
            end
    
            self.xOut = s(zOut);
        end
        
        function self = SetupEigenvalueProblem(self) 
            % The eigenvalue problem will be solved using N2 and N2z, so
            % now we need transformations to project them onto the
            % stretched grid
            T_zCheb_sLobatto = InternalModesSpectral.ChebyshevTransformForGrid(self.zLobatto, self.z_xLobatto);
            self.N2_xLobatto = T_zCheb_sLobatto(self.N2_zCheb);
            self.N2z_xLobatto = T_zCheb_sLobatto(self.Diff1_zCheb(self.N2_zCheb));
            
            Ls = max(self.xLobatto)-min(self.xLobatto);
            self.Diff1_xCheb = @(v) (2/Ls)*InternalModesSpectral.DifferentiateChebyshevVector( v );
            [self.T_xLobatto,self.Tx_xLobatto,self.Txx_xLobatto] = InternalModesSpectral.ChebyshevPolynomialsOnGrid( self.xLobatto, length(self.xLobatto) );
            [self.T_xCheb_zOut, ~] = InternalModesSpectral.ChebyshevTransformForGrid(self.xLobatto, self.xOut);
            
            % We use that \int_{-1}^1 T_n(x) dx = \frac{(-1)^n + 1}{1-n^2}
            % for all n, except n=1, where the integral is zero.
            np = (0:(self.nEVP-1))';
            self.Int_xCheb = -(1+(-1).^np)./(np.*np-1);
            self.Int_xCheb(2) = 0;
            self.Int_xCheb = Ls/2*self.Int_xCheb;
            
%             fprintf(' The eigenvalue problem will be solved with %d points.\n', length(self.xLobatto));
        end
    end
    
    methods (Access = private)
        function [F,G,h] = ModesFromGEPDensitySpectral(self,A,B)
            % This function is an intermediary used by ModesAtFrequency and
            % ModesAtWavenumber to establish the various norm functions.
            hFromLambda = @(lambda) 1.0 ./ lambda;
            GOutFromGCheb = @(G_cheb,h) self.T_xCheb_zOut(G_cheb);
            FOutFromGCheb = @(G_cheb,h) h * self.N2 .* self.T_xCheb_zOut(self.Diff1_xCheb(G_cheb));
            GFromGCheb = @(G_cheb,h) InternalModesSpectral.ifct(G_cheb);
            FFromGCheb = @(G_cheb,h) h * self.N2_xLobatto .* InternalModesSpectral.ifct( self.Diff1_xCheb(G_cheb) );
            GNorm = @(Gj) abs(Gj(1)*Gj(1) + sum(self.Int_xCheb .*InternalModesSpectral.fct((1/self.g) * (1 - self.f0*self.f0./self.N2_xLobatto) .* Gj .^ 2)));
            FNorm = @(Fj) abs(sum(self.Int_xCheb .*InternalModesSpectral.fct((1/self.Lz) * (Fj.^ 2)./self.N2_xLobatto)));
            [F,G,h] = ModesFromGEP(self,A,B,hFromLambda,GFromGCheb,FFromGCheb,GNorm,FNorm,GOutFromGCheb,FOutFromGCheb);
        end
    end
    
end
