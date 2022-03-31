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
        function [A,B] = EigenmatricesForWavenumber(self, k )
            T = self.T_xLobatto;
            Tz = self.Tx_xLobatto;
            Tzz = self.Txx_xLobatto;
                         
            A = diag(self.N2_xLobatto .* self.N2_xLobatto)*Tzz + diag(self.N2z_xLobatto)*Tz - k*k*T;
            B = diag( (self.f0*self.f0 - self.N2_xLobatto)/self.g )*T;
            
            [A,B] = self.ApplyBoundaryConditions(A,B);
        end

        function [A,B] = EigenmatricesForFrequency(self, omega )
            T = self.T_xLobatto;
            Tz = self.Tx_xLobatto;
            Tzz = self.Txx_xLobatto;
                         
            A = diag(self.N2_xLobatto .* self.N2_xLobatto)*Tzz + diag(self.N2z_xLobatto)*Tz;
            B = diag( (omega*omega - self.N2_xLobatto)/self.g )*T;
            
            [A,B] = self.ApplyBoundaryConditions(A,B);
        end

        function [A,B] = ApplyBoundaryConditions(self,A,B)
            T = self.T_xLobatto;
            Tz = self.Tx_xLobatto;
            n = self.nEVP;
            
            iSurface = 1;
            iBottom = n;
            
            switch self.lowerBoundary
                case LowerBoundary.freeSlip
                    A(iBottom,:) = T(iBottom,:);
                    B(iBottom,:) = 0;
                case LowerBoundary.noSlip
                    A(iBottom,:) = Tz(iBottom,:);
                    B(iBottom,:) = 0;
                case LowerBoundary.none
                otherwise
                    error('Unknown boundary condition');
            end
            
            % G=0 or N^2 G_s = \frac{1}{h_j} G at the surface, depending on the BC
            switch self.upperBoundary
                case UpperBoundary.freeSurface
                    % G_z = \frac{1}{h_j} G at the surface
                    A(iSurface,:) = self.N2_xLobatto(iSurface) * Tz(iSurface,:);
                    B(iSurface,:) = T(iSurface,:);
                case UpperBoundary.rigidLid
                    A(iSurface,:) = T(iSurface,:);
                    B(iSurface,:) = 0;
                case UpperBoundary.none
                otherwise
                    error('Unknown boundary condition');
            end
        end
        
    end
    
    methods (Access = protected)
        
        function self = InitializeWithGrid(self, rho, zIn)
            InitializeWithGrid@InternalModesSpectral(self,rho,zIn);
            
            % This is our new (s)treched grid...we need to make s(z) be
            % monotonic!
            s = @(z) (-self.g/self.rho0)*self.rho_function(z) + self.g;
            
            self.xDomain = [s(self.zMin) s(self.zMax)];
            
            % if the data is noisy, let's make sure that there are no
            % variations below the grid scale. Specifically, each grid
            % point should uniquely map to a z-grid point.
            if any(diff(rho)./diff(zIn) > 0)
                n = self.nEVP;
                xLobatto = ((self.xMax-self.xMin)/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + self.xMin;
                xIn = s(zIn);
                y = discretize(xLobatto,xIn);
                while (length(unique(y)) < length(xLobatto))
                    n = n-1;
                    xLobatto = ((self.xMax-self.xMin)/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + self.xMin;
                    y = discretize(xLobatto,xIn);
                end
                
                K = 6; % cubic spline
                z_knot = InterpolatingSpline.KnotPointsForPoints([zIn(1);zIn(unique(y)+1)],K,1);
                rho_interpolant = ConstrainedSpline(zIn,rho,K,z_knot,NormalDistribution(1),struct('global',ShapeConstraint.monotonicDecreasing));
                
                self.rho_function = rho_interpolant;
                self.N2_function = (-self.g/self.rho0)*diff(self.rho_function);
                
                if self.shouldShowDiagnostics == 1
                    fprintf('Data was found to be noise. Creating a %d-order monotonic smoothing spline using %d knot points.\n', K, n);
                end
            end
        end
        
        function self = SetupEigenvalueProblem(self)     
            % Create a stretched grid from the density function
            self.x_function = @(z) (-self.g/self.rho0)*self.rho_function(z) + self.g;
            
            N2z_function = diff(self.N2_function);
            self.N2z_xLobatto = N2z_function(self.z_xLobatto);

            self.hFromLambda = @(lambda) 1.0 ./ lambda;
            self.GOutFromGCheb = @(G_cheb,h) self.T_xCheb_zOut(G_cheb);
            self.FOutFromGCheb = @(G_cheb,h) h * self.N2 .* self.T_xCheb_zOut(self.Diff1_xCheb(G_cheb));
            self.GFromGCheb = @(G_cheb,h) InternalModesSpectral.ifct(G_cheb);
            self.FFromGCheb = @(G_cheb,h) h * self.N2_xLobatto .* InternalModesSpectral.ifct( self.Diff1_xCheb(G_cheb) );
            self.GNorm = @(Gj) abs(Gj(1)*Gj(1) + sum(self.Int_xCheb .*InternalModesSpectral.fct((1/self.g) * (1 - self.f0*self.f0./self.N2_xLobatto) .* Gj .^ 2)));
            self.FNorm = @(Fj) abs(sum(self.Int_xCheb .*InternalModesSpectral.fct((1/self.Lz) * (Fj.^ 2)./self.N2_xLobatto)));
        end
    end
    
end
