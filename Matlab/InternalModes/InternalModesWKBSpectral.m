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
                        
            A = diag(self.N2_xLobatto)*Tzz + diag(self.Nz_xLobatto)*Tz - k*k*T;
            B = diag( (self.f0*self.f0 - self.N2_xLobatto)/self.g )*T;
            
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
                    % N*G_s = \frac{1}{h_j} G at the surface
                    A(1,:) = sqrt(self.N2_xLobatto(1)) * Tz(1,:);
                    B(1,:) = T(1,:);
                case UpperBoundary.rigidLid
                    A(1,:) = T(1,:);
                    B(1,:) = 0;
                case UpperBoundary.none
                otherwise
                    error('Unknown boundary condition');
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
                    % N*G_s = \frac{1}{h_j} G at the surface
                    A(1,:) = sqrt(self.N2_xLobatto(1)) * Tz(1,:);
                    B(1,:) = T(1,:);
                case UpperBoundary.rigidLid
                    A(1,:) = T(1,:);
                    B(1,:) = 0;
                case UpperBoundary.none
                otherwise
                    error('Unknown boundary condition');
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
        
        function self = InitializeWithGrid(self, rho, zIn)
            InitializeWithGrid@InternalModesSpectral(self,rho,zIn);
            rho=flip(rho); zIn=flip(zIn);
            
            N_function = sqrt(self.N2_function);
            s = cumsum(N_function);
            self.xDomain = [s(self.zMin) s(self.zMax)];
            
            n = self.nEVP;
            
            self.xLobatto = ((self.xMax-self.xMin)/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + self.xMin;
            [self.z_xLobatto, self.xOut] = InternalModesSpectral.StretchedGridFromCoordinate( s, self.xLobatto, self.zDomain, self.z);
            y = discretize(self.z_xLobatto,zIn);
            
                        
            % The goal here is to eliminate variations below the WKB grid
            % scale. So, we don't allow knot points 
            while (length(unique(y)) < length(self.z_xLobatto))
                n = (n-1)/2+1;
                self.xLobatto = ((self.xMax-self.xMin)/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + self.xMin;
                [self.z_xLobatto, self.xOut] = InternalModesSpectral.StretchedGridFromCoordinate( s, self.xLobatto, self.zDomain, self.z);
                y = discretize(self.z_xLobatto,zIn);
            end
            

            K = 4; % cubic spline

            
            
            
%             z_knotIn = InterpolatingSpline.KnotPointsForPoints(zIn,K,1);
            z_knot = InterpolatingSpline.KnotPointsForPoints([zIn(1);zIn(unique(y)+1)],K,1);
            z_knot = InterpolatingSpline.KnotPointsForPoints(self.z_xLobatto,K,1);
            rho_interpolant = ConstrainedSpline(zIn,rho,K,z_knot,NormalDistribution(1),struct('global',ShapeConstraint.monotonicDecreasing));
            
            self.rho_function = rho_interpolant;
            self.N2_function = (-self.g/self.rho0)*diff(self.rho_function);
        end
        
        function self = SetupEigenvalueProblem(self)       
            % Create the stretched WKB grid
            N_function = sqrt(self.N2_function);
            s = cumsum(N_function);
            self.xDomain = [s(self.zMin) s(self.zMax)];
            self.xLobatto = ((self.xMax-self.xMin)/2)*( cos(((0:self.nEVP-1)')*pi/(self.nEVP-1)) + 1) + self.xMin;

            [self.z_xLobatto, self.xOut] = InternalModesSpectral.StretchedGridFromCoordinate( s, self.xLobatto, self.zDomain, self.z);

            
            % The eigenvalue problem will be solved using N2 and N2z, so
            % now we need transformations to project them onto the
            % stretched grid
            self.N2_xLobatto = self.N2_function(self.z_xLobatto);
            Nz_function = diff(N_function);
            self.Nz_xLobatto = Nz_function(self.z_xLobatto);
            
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
            
            if self.shouldShowDiagnostics == 1
                fprintf(' The eigenvalue problem will be solved with %d points.\n', length(self.xLobatto));
            end
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
