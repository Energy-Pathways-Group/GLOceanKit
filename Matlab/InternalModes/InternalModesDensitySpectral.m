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
        function self = SetupEigenvalueProblem(self)
            
            % Check if we an even create a density coordinate system
            [flag, dTotalVariation, rho_zCheb_new, rho_zLobatto_new, rhoz_zCheb, rhoz_zLobatto] = InternalModesSpectral.CheckIfReasonablyMonotonic(self.zLobatto, self.rho_zCheb, self.rho_zLobatto, -(self.rho0/self.g)*self.N2_zCheb, -(self.rho0/self.g)*self.N2_zLobatto);
            if flag == 1
                self.rho_zCheb = rho_zCheb_new;
                self.rho_zLobatto = rho_zLobatto_new;
                self.N2_zCheb= -(self.g/self.rho0)*rhoz_zCheb;
                self.N2_zLobatto = -(self.g/self.rho0)*rhoz_zLobatto;
                
                self.rho = self.T_zCheb_zOut(self.rho_zCheb);
                self.N2 = self.T_zCheb_zOut(self.N2_zCheb);
                
                fprintf('The density function was not monotonically decreasing and zeroing out overturns resulted in a change in total variation of %.2g percent. We used this new density function for the computation and will proceed.\n', dTotalVariation*100);
            elseif flag == 2
                error('The density function was not monotonically decreasing and zeroing out overturns resulted in a change in total variation of %.2g percent. We are unable to create a WKB stretched coordinate system.\n', dTotalVariation*100);
            end
            
            
            % Create a stretched grid from the density function
            s = @(z) -self.g*self.rho_function(z)/self.rho0 + self.g;
            
            n = self.nEVP;
            range = [s(max(self.zLobatto)), s(min(self.zLobatto))];
            sMin = min(range);
            sMax = max(range);
            Ls = sMax-sMin;
            self.xLobatto = (Ls/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + sMin;
            [self.z_xLobatto, self.xOut] = InternalModesSpectral.StretchedGridFromCoordinate( s, self.xLobatto, self.zLobatto, self.z);
            
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
