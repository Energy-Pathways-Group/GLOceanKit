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
            varargin{end+1} = 'requiresMonotonicDensity';
            varargin{end+1} = 1;
            self@InternalModesSpectral(rho,z_in,z_out,latitude, varargin{:});
        end
                    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computation of the modes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [F,G,h,omega,F2,N2G2] = ModesAtWavenumber(self, k )
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
            
            if nargout == 6
                [F,G,h,F2,N2G2] = self.ModesFromGEPDensitySpectral(A,B);
            else
                [F,G,h] = self.ModesFromGEPDensitySpectral(A,B);
            end
            omega = self.omegaFromK(h,k);
        end
        
        function [F,G,h,k,F2,N2G2] = ModesAtFrequency(self, omega )
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
            
            if nargout == 6
                [F,G,h,F2,N2G2] = self.ModesFromGEPDensitySpectral(A,B);
            else
                [F,G,h] = self.ModesFromGEPDensitySpectral(A,B);
            end
            k = self.kFromOmega(h,omega);
        end

    end
    
    methods (Access = protected)        
        function self = SetupEigenvalueProblem(self)     
            % Create a stretched grid from the density function
            s = @(z) (-self.g/self.rho0)*self.rho_function(z) + self.g;
            self.xDomain = [s(self.zMin) s(self.zMax)];
            self.xLobatto = ((self.xMax-self.xMin)/2)*( cos(((0:self.nEVP-1)')*pi/(self.nEVP-1)) + 1) + self.xMin;
            
            [self.z_xLobatto, self.xOut] = InternalModesSpectral.StretchedGridFromCoordinate( s, self.xLobatto, self.zDomain, self.z);
            
            % The eigenvalue problem will be solved using N2 and N2z, so
            % now we need transformations to project them onto the
            % stretched grid
            self.N2_xLobatto = self.N2_function(self.z_xLobatto);
            N2z_function = diff(self.N2_function);
            self.N2z_xLobatto = N2z_function(self.z_xLobatto);
            
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
            
            if self.shouldShowDiagnostics == 1
                fprintf(' The eigenvalue problem will be solved with %d points.\n', length(self.xLobatto));
            end
        end
    end
    
    methods (Access = private)
        function [F,G,h,F2,N2G2] = ModesFromGEPDensitySpectral(self,A,B)
            % This function is an intermediary used by ModesAtFrequency and
            % ModesAtWavenumber to establish the various norm functions.
            hFromLambda = @(lambda) 1.0 ./ lambda;
            GOutFromGCheb = @(G_cheb,h) self.T_xCheb_zOut(G_cheb);
            FOutFromGCheb = @(G_cheb,h) h * self.N2 .* self.T_xCheb_zOut(self.Diff1_xCheb(G_cheb));
            GFromGCheb = @(G_cheb,h) InternalModesSpectral.ifct(G_cheb);
            FFromGCheb = @(G_cheb,h) h * self.N2_xLobatto .* InternalModesSpectral.ifct( self.Diff1_xCheb(G_cheb) );
            GNorm = @(Gj) abs(Gj(1)*Gj(1) + sum(self.Int_xCheb .*InternalModesSpectral.fct((1/self.g) * (1 - self.f0*self.f0./self.N2_xLobatto) .* Gj .^ 2)));
            FNorm = @(Fj) abs(sum(self.Int_xCheb .*InternalModesSpectral.fct((1/self.Lz) * (Fj.^ 2)./self.N2_xLobatto)));
            if nargout == 5
                [F,G,h,F2,N2G2] = ModesFromGEP(self,A,B,hFromLambda,GFromGCheb,FFromGCheb,GNorm,FNorm,GOutFromGCheb,FOutFromGCheb);
            else
                [F,G,h] = ModesFromGEP(self,A,B,hFromLambda,GFromGCheb,FFromGCheb,GNorm,FNorm,GOutFromGCheb,FOutFromGCheb);
            end
        end
    end
    
end
