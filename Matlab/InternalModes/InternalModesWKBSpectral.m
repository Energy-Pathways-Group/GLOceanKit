classdef InternalModesWKBSpectral < InternalModesSpectral
    % This class solves the vertical eigenvalue problem on a WKB stretched
    % density coordinate grid using Chebyshev polynomials.
    %
    % See InternalModesBase for basic usage information.
    %
    % This class uses the coordinate
    %   s = \int_{-Lz}^0 \sqrt(-(g/rho0)*rho_z) dz
    % to solve the EVP.
    %
    % Internally, xLobatto is the stretched WKB coordinate on a
    % Chebyshev extrema/Lobatto grid. This is the grid upon which the
    % eigenvalue problem is solved, and therefore the class uses the
    % superclass properties denoted with 'x'.
    %
    %   See also INTERNALMODES, INTERNALMODESBASE, INTERNALMODESSPECTRAL,
    %   INTERNALMODESDENSITYSPECTRAL, and INTERNALMODESFINITEDIFFERENCE.
    %
    %   Jeffrey J. Early
    %   jeffrey@jeffreyearly.com
    %
    %   March 14th, 2017        Version 1.0
    
    properties %(Access = private)            
        Nz_xLobatto     	% (d/dz)N on the xiLobatto grid   
    end
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = InternalModesWKBSpectral(rho, z_in, z_out, latitude, varargin)
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
                        
            A = diag(self.N2_xLobatto)*Tzz + diag(self.Nz_xLobatto)*Tz - k*k*T;
            B = diag( (self.f0*self.f0 - self.N2_xLobatto)/self.g )*T;
            
            % Lower boundary is rigid, G=0
            A(n,:) = T(n,:);
            B(n,:) = 0;
            
            % G=0 or N*G_s = \frac{1}{h_j} G at the surface, depending on the BC
            if self.upperBoundary == UpperBoundary.freeSurface
                A(1,:) = sqrt(self.N2_xLobatto(1)) * Tz(1,:);
                B(1,:) = T(1,:);
            elseif self.upperBoundary == UpperBoundary.rigidLid
                A(1,:) = T(1,:);
                B(1,:) = 0;
            end
            
            if nargout == 6
                [F,G,h,F2,N2G2] = self.ModesFromGEPWKBSpectral(A,B);
            else
                [F,G,h] = self.ModesFromGEPWKBSpectral(A,B);
            end
            omega = self.omegaFromK(h,k);
        end
        
        function [F,G,h,k,F2,N2G2] = ModesAtFrequency(self, omega )
            self.gridFrequency = omega;
            
            T = self.T_xLobatto;
            Tz = self.Tx_xLobatto;
            Tzz = self.Txx_xLobatto;
            n = self.nEVP;
            
            A = diag(self.N2_xLobatto)*Tzz + diag(self.Nz_xLobatto)*Tz;
            B = diag( (omega*omega - self.N2_xLobatto)/self.g )*T;
            
            % Lower boundary is rigid, G=0
            A(n,:) = T(n,:);
            B(n,:) = 0;
            
            % G=0 or N*G_s = \frac{1}{h_j} G at the surface, depending on the BC
            if self.upperBoundary == UpperBoundary.freeSurface
                A(1,:) = sqrt(self.N2_xLobatto(1)) * Tz(1,:);
                B(1,:) = T(1,:);
            elseif self.upperBoundary == UpperBoundary.rigidLid
                A(1,:) = T(1,:);
                B(1,:) = 0;
            end
            
            if nargout == 6
                [F,G,h,F2,N2G2] = self.ModesFromGEPWKBSpectral(A,B);
            else
                [F,G,h] = self.ModesFromGEPWKBSpectral(A,B);
            end
            k = self.kFromOmega(h,omega);
        end
 
    end
    
    methods (Access = protected)        
        function self = SetupEigenvalueProblem(self)
            
            % Check if we an even create a WKB coordinate system
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
            
            % Create the stretched WKB grid
            N_zLobatto = sqrt(self.N2_zLobatto);
            N_zCheb = InternalModesSpectral.fct(N_zLobatto);  
            x_zCheb = (self.Lz/2)*InternalModesSpectral.IntegrateChebyshevVector(N_zCheb);
            
            s = @(z) InternalModesSpectral.ValueOfFunctionAtPointOnGrid(z,self.zLobatto,x_zCheb);
            self.xDomain = [s(self.zMin) s(self.zMax)];
            self.xLobatto = ((self.xMax-self.xMin)/2)*( cos(((0:self.nEVP-1)')*pi/(self.nEVP-1)) + 1) + self.xMin;

            % We need to be able to create a reasonable stretched grid...
            % if we can't, we will throw an exception
            try
                [self.z_xLobatto, self.xOut] = InternalModesSpectral.StretchedGridFromCoordinate( s, self.xLobatto, self.zLobatto, self.z);
            catch ME
                switch ME.identifier
                    case 'MATLAB:griddedInterpolant:NonUniqueCompVecsPtsErrId'
                        causeException = MException('StretchedGridFromCoordinate:NonMonotonicFunction','The density function must be strictly monotonically decreasing in order to create a unqiue wkb grid. Unable to proceed. You consider trying InternalModesSpectral, which has no such restriction because it uses the z-coordinate.');
                        ME = addCause(ME,causeException);
                end
                rethrow(ME)
            end
                        
            % The eigenvalue problem will be solved using N2 and N2z, so
            % now we need transformations to project them onto the
            % stretched grid
            T_zCheb_xiLobatto = InternalModesSpectral.ChebyshevTransformForGrid(self.zLobatto, self.z_xLobatto);
            self.N2_xLobatto = T_zCheb_xiLobatto(self.N2_zCheb);
            self.Nz_xLobatto = T_zCheb_xiLobatto(self.Diff1_zCheb(N_zCheb));
            
            Lxi = max(self.xLobatto) - min(self.xLobatto);
            self.Diff1_xCheb = @(v) (2/Lxi)*InternalModesSpectral.DifferentiateChebyshevVector( v );
            [self.T_xLobatto,self.Tx_xLobatto,self.Txx_xLobatto] = InternalModesSpectral.ChebyshevPolynomialsOnGrid( self.xLobatto, length(self.xLobatto) );
            
            self.T_xCheb_zOut = InternalModesSpectral.ChebyshevTransformForGrid(self.xLobatto, self.xOut);
            
            % We use that \int_{-1}^1 T_n(x) dx = \frac{(-1)^n + 1}{1-n^2}
            % for all n, except n=1, where the integral is zero.
            np = (0:(self.nEVP-1))';
            self.Int_xCheb = -(1+(-1).^np)./(np.*np-1);
            self.Int_xCheb(2) = 0;
            self.Int_xCheb = Lxi/2*self.Int_xCheb;
            
%             fprintf(' The eigenvalue problem will be solved with %d points.\n', length(self.xLobatto));
        end
        
        
        
    end
    
    methods (Access = private)             
        function [F,G,h,F2,N2G2] = ModesFromGEPWKBSpectral(self,A,B)
            % This function is an intermediary used by ModesAtFrequency and
            % ModesAtWavenumber to establish the various norm functions.
            hFromLambda = @(lambda) 1.0 ./ lambda;
            GOutFromGCheb = @(G_cheb,h) self.T_xCheb_zOut(G_cheb);
            FOutFromGCheb = @(G_cheb,h) h * sqrt(self.N2) .* self.T_xCheb_zOut(self.Diff1_xCheb(G_cheb));
            GFromGCheb = @(G_cheb,h) InternalModesSpectral.ifct(G_cheb);
            FFromGCheb = @(G_cheb,h) h * sqrt(self.N2_xLobatto) .* InternalModesSpectral.ifct( self.Diff1_xCheb(G_cheb) );
            GNorm = @(Gj) abs(Gj(1)*Gj(1) + sum(self.Int_xCheb .*InternalModesSpectral.fct((1/self.g) * (self.N2_xLobatto - self.f0*self.f0) .* ( self.N2_xLobatto.^(-0.5) ) .* Gj .^ 2)));
            FNorm = @(Fj) abs(sum(self.Int_xCheb .*InternalModesSpectral.fct((1/self.Lz) * (Fj.^ 2) .* ( self.N2_xLobatto.^(-0.5) ))));
            if nargout == 5
                [F,G,h,F2,N2G2] = ModesFromGEP(self,A,B,hFromLambda,GFromGCheb,FFromGCheb,GNorm,FNorm,GOutFromGCheb,FOutFromGCheb);
            else
                [F,G,h] = ModesFromGEP(self,A,B,hFromLambda,GFromGCheb,FFromGCheb,GNorm,FNorm,GOutFromGCheb,FOutFromGCheb);
            end
        end
    end
    
end
