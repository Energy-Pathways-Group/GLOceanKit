classdef InternalModesAdaptiveSpectral < InternalModesSpectral
    % This class solves the vertical eigenvalue problem on a WKB stretched
    % density coordinate grid using Chebyshev polynomials.
    %
    % See InternalModesBase for basic usage information.
    %
    % This class uses the coordinate
    %   s = \int_{-Lz}^0 \sqrt(-(g/rho0)*rho_z) dz
    % to solve the EVP.
    %
    % Internally, sLobatto is the stretched WKB coordinate on a
    % Chebyshev extrema/Lobatto grid. This is the grid upon which the
    % eigenvalue problem is solved, and therefore the class uses the
    % superclass properties denoted with 'x' instead of 's' when setting up
    % the eigenvalue problem.
    %
    %   See also INTERNALMODES, INTERNALMODESBASE, INTERNALMODESSPECTRAL,
    %   INTERNALMODESDENSITYSPECTRAL, and INTERNALMODESFINITEDIFFERENCE.
    %
    %   Jeffrey J. Early
    %   jeffrey@jeffreyearly.com
    %
    %   March 14th, 2017        Version 1.0
    
    properties %(Access = private)
        gridFrequency        % The value of omega used for s = int sqrt(N^2 - omega^2) dz
        xiLobatto            % stretched density coordinate, on Chebyshev extrema/Lobatto grid
        z_xiLobatto          % The value of z, at the sLobatto points
        xiOut                % desired locations of the output in s-coordinate (deduced from z_out)
        
        zBoundaries         % z-location of the boundaries (end points plus turning points).
        xiBoundaries        % xi-location of the boundaries (end points plus turning points).
        boundaryIndices     % indices of the boundaries into xiLobatto
        
        SqrtN2Omega2_xLobatto
        N2Omega2_xLobatto
    end
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = InternalModesAdaptiveSpectral(rho, z_in, z_out, latitude, varargin)
            self@InternalModesSpectral(rho,z_in,z_out,latitude, varargin{:});
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computation of the modes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [F,G,h] = ModesAtWavenumber(self, k )
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
            if strcmp(self.upperBoundary, 'free_surface')
                A(1,:) = sqrt(self.N2_xLobatto(1)) * Tz(1,:);
                B(1,:) = T(1,:);
            elseif strcmp(self.upperBoundary, 'rigid_lid')
                A(1,:) = T(1,:);
                B(1,:) = 0;
            end
            
            [F,G,h] = self.ModesFromGEPWKBSpectral(A,B);
        end
        
        function [F,G,h] = ModesAtFrequency(self, omega )
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
            if strcmp(self.upperBoundary, 'free_surface')
                A(1,:) = sqrt(self.N2_xLobatto(1)) * Tz(1,:);
                B(1,:) = T(1,:);
            elseif strcmp(self.upperBoundary, 'rigid_lid')
                A(1,:) = T(1,:);
                B(1,:) = 0;
            end
            
            [F,G,h] = self.ModesFromGEPWKBSpectral(A,B);
        end
 
    end
    
    methods (Access = protected)
        function self = InitializeWithGrid(self, rho, z_in)
            % Superclass calls this method upon initialization when it
            % determines that the input is given in gridded form.
            %
            % The superclass will initialize zLobatto and rho_lobatto; this
            % class must initialize the sLobatto, z_sLobatto and sOut.
            InitializeWithGrid@InternalModesSpectral(self, rho, z_in);
            
            self.InitializeWKBGridWithFrequency(0)
        end
        

        function self = InitializeWithFunction(self, rho, zMin, zMax, zOut)
            % Superclass calls this method upon initialization when it
            % determines that the input is given in functional form.
            %
            % The superclass will initialize zLobatto and rho_lobatto; this
            % class must initialize the sLobatto, z_sLobatto and sOut.
            InitializeWithFunction@InternalModesSpectral(self, rho, zMin, zMax, zOut);
            
            self.InitializeWKBGridWithFrequency(3*self.f0)
        end
        
        function InitializeWKBGridWithFrequency(self,omega)
            self.gridFrequency = omega;
            
            % Create the stretched grid \xi
            N2Omega2_zLobatto = InternalModesSpectral.ifct(self.N2_zCheb) - self.gridFrequency*self.gridFrequency;
            xi_zLobatto = cumtrapz(self.zLobatto,sqrt(abs(N2Omega2_zLobatto)));
            
            % Now find (roughly) the turning points, if any
            a = N2Omega2_zLobatto; a(a>=0) = 1; a(a<0) = 0;
            turningIndices = find(diff(a)~=0);
            nTP = length(turningIndices);
            zTP = zeros(nTP,1);
            for i=1:nTP
                fun = @(z) interp1(self.zLobatto,N2Omega2_zLobatto,z,'spline');
                zTP(i) = fzero(fun,self.zLobatto(turningIndices(i)));
            end  
            self.zBoundaries = [self.zLobatto(1); zTP; self.zLobatto(end)];
            self.xiBoundaries = interp1(self.zLobatto,xi_zLobatto,self.zBoundaries);
            
            nBoundaries = length(self.xiBoundaries);
            nEquations = nTP+1;
            
            % We will be coupling nTP+1 EVPs together. We need to
            % distribute the user requested points to each of these EVPs.
            % For this first draft, we simply evenly distribute the points.
            nPoints = floor((self.nEVP-nBoundaries)/nEquations);
            nInteriorPoints = nPoints*ones(nEquations,1);
            nInteriorPoints(end) = nInteriorPoints(end) + (self.nEVP-nBoundaries) - nPoints*nEquations; % add the extra points to the end
            self.boundaryIndices = cumsum( [1; nInteriorPoints+1] );
            
            self.xiLobatto = zeros(self.nEVP,1);
            for i=1:nEquations
                Lxi = max(self.xiBoundaries(i+1),self.xiBoundaries(i)) - min(self.xiBoundaries(i+1),self.xiBoundaries(i));
                n = nInteriorPoints(i)+2;
                self.xiLobatto(self.boundaryIndices(i):self.boundaryIndices(i+1)) = (Lxi/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + min(self.xiBoundaries(i+1),self.xiBoundaries(i));
            end
            
            % Now we need z on the \xi grid
            self.z_xiLobatto = interp1(xi_zLobatto, self.zLobatto, self.xiLobatto, 'spline');
            
            % and z_out on the \xi grid
            self.xiOut = interp1(self.zLobatto, xi_zLobatto, self.z, 'spline');
        end
        
        function self = SetupEigenvalueProblem(self)
            nEquations = length(self.zBoundaries)-1;
            
            % We will use the stretched grid to solve the eigenvalue
            % problem.
            self.xLobatto = self.xiLobatto;
            self.SqrtN2Omega2_xLobatto = zeros(size(self.xLobatto));
            self.N2Omega2_xLobatto = zeros(size(self.xLobatto));
            
            for i=1:nEquations
                L = abs(self.zBoundaries(i)-self.zBoundaries(i+1));
                zMin = min(self.zBoundaries(i),self.zBoundaries(i+1));
                n = self.nGrid;
                zEqLobatto = (L/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + zMin;
                rho_zEqLobatto = self.rho_function(zEqLobatto);
                rho_zEqCheb = InternalModesSpectral.fct(rho_zEqLobatto);
                rho_zEqCheb = self.SetNoiseFloorToZero(rho_zEqCheb);
                N2_zEqCheb= -(self.g/self.rho0)*(2/L)*InternalModesSpectral.DifferentiateChebyshevVector(rho_zEqCheb);
                N2_zEqLobatto = InternalModesSpectral.ifct(N2_zEqCheb);
                SqrtN2Omega2_zEqLobatto = sqrt(abs(N2_zEqLobatto - self.gridFrequency*self.gridFrequency));
                SqrtN2Omega2_zEqCheb = InternalModesSpectral.fct(SqrtN2Omega2_zEqLobatto);
                N2Omega2_zEqCheb = N2_zEqCheb;
                N2Omega2_zEqCheb(1) = N2Omega2_zEqCheb(1)-self.gridFrequency*self.gridFrequency;
                
                indices = self.boundaryIndices(i):self.boundaryIndices(i+1);
                zEq_xiLobatto = self.z_xiLobatto( indices );
                T_zEqCheb_xiLobatto = InternalModesSpectral.ChebyshevTransformForGrid(zEqLobatto, zEq_xiLobatto);
                
                self.SqrtN2Omega2_xLobatto(indices) = T_zEqCheb_xiLobatto(SqrtN2Omega2_zEqCheb);
                self.N2Omega2_xLobatto(indices) = T_zEqCheb_xiLobatto(N2Omega2_zEqCheb);
            end
            
            %%%% Current stopping point, Sept 7, 2017

            
            % The eigenvalue problem will be solved using N2 and N2z, so
            % now we need transformations to project them onto the
            % stretched grid
            T_zCheb_xiLobatto = InternalModesSpectral.ChebyshevTransformForGrid(self.zLobatto, self.z_xiLobatto);
            self.N2_xLobatto = T_zCheb_xiLobatto(self.N2_zCheb);
            self.Nz_xLobatto = T_zCheb_xiLobatto(self.Diff1_zCheb(self.N_zCheb));
            
            Lxi = max(self.xiLobatto) - min(self.xiLobatto);
            self.Diff1_xCheb = @(v) (2/Lxi)*InternalModesSpectral.DifferentiateChebyshevVector( v );
            [self.T_xLobatto,self.Tx_xLobatto,self.Txx_xLobatto] = InternalModesSpectral.ChebyshevPolynomialsOnGrid( self.xiLobatto, length(self.xiLobatto) );
            
            self.T_xCheb_zOut = InternalModesSpectral.ChebyshevTransformForGrid(self.xiLobatto, self.xiOut);
            
            % We use that \int_{-1}^1 T_n(x) dx = \frac{(-1)^n + 1}{1-n^2}
            % for all n, except n=1, where the integral is zero.
            np = (0:(self.nEVP-1))';
            self.Int_xCheb = -(1+(-1).^np)./(np.*np-1);
            self.Int_xCheb(2) = 0;
            self.Int_xCheb = Lxi/2*self.Int_xCheb;
        end
    end
    
    methods (Access = private)             
        function [F,G,h] = ModesFromGEPWKBSpectral(self,A,B)
            % This function is an intermediary used by ModesAtFrequency and
            % ModesAtWavenumber to establish the various norm functions.
            hFromLambda = @(lambda) 1.0 ./ lambda;
            GOutFromGCheb = @(G_cheb,h) self.T_xCheb_zOut(G_cheb);
            FOutFromGCheb = @(G_cheb,h) h * sqrt(self.N2) .* self.T_xCheb_zOut(self.Diff1_xCheb(G_cheb));
            GFromGCheb = @(G_cheb,h) InternalModesSpectral.ifct(G_cheb);
            FFromGCheb = @(G_cheb,h) h * sqrt(self.N2_xLobatto) .* InternalModesSpectral.ifct( self.Diff1_xCheb(G_cheb) );
            GNorm = @(Gj) abs(sum(self.Int_xCheb .*InternalModesSpectral.fct((1/self.g) * (self.N2_xLobatto - self.f0*self.f0) .* ( self.N2_xLobatto.^(-0.5) ) .* Gj .^ 2)));
            FNorm = @(Fj) abs(sum(self.Int_xCheb .*InternalModesSpectral.fct((1/self.Lz) * (Fj.^ 2) .* ( self.N2_xLobatto.^(-0.5) ))));
            [F,G,h] = ModesFromGEP(self,A,B,hFromLambda,GFromGCheb,FFromGCheb,GNorm,FNorm,GOutFromGCheb,FOutFromGCheb);
        end
    end
    
end
