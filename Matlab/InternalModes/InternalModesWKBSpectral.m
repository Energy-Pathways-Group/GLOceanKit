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
        Nz_function
        Nz_xLobatto     	% (d/dz)N on the xiLobatto grid   
        requiresMonotonicDensitySetting = 1
    end
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = InternalModesWKBSpectral(options)
            arguments
                options.rho = ''
                options.N2 function_handle = @disp
                options.zIn (:,1) double = []
                options.zOut (:,1) double = []
                options.latitude (1,1) double = 33
                options.rho0 (1,1) double {mustBePositive} = 1025
                options.nModes (1,1) double = 0
                options.nEVP = 512;
                options.rotationRate (1,1) double = 7.2921e-5;
                options.g (1,1) double = 9.81
            end
            self@InternalModesSpectral(rho=options.rho,N2=options.N2,zIn=options.zIn,zOut=options.zOut,latitude=options.latitude,rho0=options.rho0,nModes=options.nModes,nEVP=options.nEVP,rotationRate=options.rotationRate,g=options.g);

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
                        
            A = diag(self.N2_xLobatto)*Tzz + diag(self.Nz_xLobatto)*Tz - k*k*T;
            B = diag( (self.f0*self.f0 - self.N2_xLobatto)/self.g )*T;
            
            [A,B] = self.ApplyBoundaryConditions(A,B);
        end
        
        function [A,B] = EigenmatricesForFrequency(self, omega )
            T = self.T_xLobatto;
            Tz = self.Tx_xLobatto;
            Tzz = self.Txx_xLobatto;
            
            A = diag(self.N2_xLobatto)*Tzz + diag(self.Nz_xLobatto)*Tz;
            B = diag( (omega*omega - self.N2_xLobatto)/self.g )*T;
            
            [A,B] = self.ApplyBoundaryConditions(A,B);
        end

        function [A,B] = EigenmatricesForMDAModes(self )
            T = self.T_xLobatto;
            Tz = self.Tx_xLobatto;
            Tzz = self.Txx_xLobatto;
            n = self.nEVP;

            A = diag(self.N2_xLobatto)*Tzz + diag(self.Nz_xLobatto)*Tz;
            B = diag( - self.N2_xLobatto/self.g )*T;

            % upper-boundary
            A(1,:) = Tz(1,:); %-Tz(n,:);
            B(1,:) = 0 ;%1/self.Lz; %0*T(n,:);
            self.upperBoundary = UpperBoundary.mda;

            % lower-boundary
            A(n,:) = Tz(n,:); %self.Lz*Tz(n,:)-T(n,:);
            B(n,:) = 0; %1/self.Lz; %0*T(n,:);
            self.lowerBoundary = LowerBoundary.mda;
        end

        function [A,B] = EigenmatricesForGeostrophicGModes(self, k )
            T = self.T_xLobatto;
            Tz = self.Tx_xLobatto;
            Tzz = self.Txx_xLobatto;
            n = self.nEVP;

            A = diag(self.N2_xLobatto)*Tzz + diag(self.Nz_xLobatto)*Tz;
            B = diag( (- self.N2_xLobatto)/self.g )*T;

            A(n,:) = T(n,:);
            B(n,:) = 0;

            if (self.g/(self.f0*self.f0))*(k*k) < 1
                A(1,:) = sqrt(self.N2_xLobatto(1))*Tz(1,:) + (self.g/(self.f0*self.f0))*(k*k)*T(1,:);
            else
                A(1,:) = sqrt(self.N2_xLobatto(1))*(self.f0*self.f0)/(self.g *k*k)* Tz(1,:) + T(1,:);
            end
            B(1,:) = 0;
        end
                
        function [A,B] = ApplyBoundaryConditions(self,A,B)
            T = self.T_xLobatto;
            Tz = self.Tx_xLobatto;
            n = self.nEVP;
            
            switch self.lowerBoundary
                case LowerBoundary.freeSlip
                    A(n,:) = T(n,:);
                    B(n,:) = 0;
                case LowerBoundary.noSlip
                    A(n,:) = Tz(n,:);
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
        end
        
    end
    
    methods (Access = protected)
        
        function self = InitializeWithGrid(self, rho, zIn)
            InitializeWithGrid@InternalModesSpectral(self,rho,zIn);

            % This is our new (s)treched grid...we need to make s(z) be
            % monotonic!
            N_function = sqrt(self.N2_function,struct('global',ShapeConstraint.positive));
            s = cumsum(N_function);
            
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
                
                K = 4; % cubic spline
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
            % Create the stretched WKB grid
            N_function = sqrt(self.N2_function);
            self.x_function = cumsum(N_function);
            
            self.Nz_function = diff(N_function);
            self.Nz_xLobatto = self.Nz_function(self.z_xLobatto);

            self.hFromLambda = @(lambda) 1.0 ./ lambda;
            self.GOutFromVCheb = @(G_cheb,h) self.T_xCheb_zOut(G_cheb);
            self.FOutFromVCheb = @(G_cheb,h) h * sqrt(self.N2) .* self.T_xCheb_zOut(self.Diff1_xCheb(G_cheb));
            self.GFromVCheb = @(G_cheb,h) InternalModesSpectral.ifct(G_cheb);
            self.FFromVCheb = @(G_cheb,h) h * sqrt(self.N2_xLobatto) .* InternalModesSpectral.ifct( self.Diff1_xCheb(G_cheb) );
            self.GNorm = @(Gj) abs(Gj(1)*Gj(1) + sum(self.Int_xCheb .*InternalModesSpectral.fct((1/self.g) * (self.N2_xLobatto - self.f0*self.f0) .* ( self.N2_xLobatto.^(-0.5) ) .* Gj .^ 2)));
            self.GeostrophicNorm = @(Gj) abs(sum(self.Int_xCheb .*InternalModesSpectral.fct((1/self.g) * self.N2_xLobatto .* ( self.N2_xLobatto.^(-0.5) ) .* Gj .^ 2)));
            self.FNorm = @(Fj) abs(sum(self.Int_xCheb .*InternalModesSpectral.fct((1/self.Lz) * (Fj.^ 2) .* ( self.N2_xLobatto.^(-0.5) ))));
        end   
    end
    
end
