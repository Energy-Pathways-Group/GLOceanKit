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
        xiLobatto            % stretched density coordinate, on Chebyshev extrema/Lobatto grid
        z_xiLobatto          % The value of z, at the sLobatto points
        xiOut                % desired locations of the output in s-coordinate (deduced from z_out)
        
        N_zCheb
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
            %
            % Superclass initializes zLobatto and rho_lobatto
            InitializeWithGrid@InternalModesSpectral(self, rho, z_in);
            
            % Create the stretched grid \xi
            N_zLobatto = sqrt(InternalModesSpectral.ifct(self.N2_zCheb));
            self.N_zCheb = InternalModesSpectral.fct(N_zLobatto);
            xi_zLobatto = cumtrapz(self.zLobatto,N_zLobatto);
            Lxi = max(xi_zLobatto)-min(xi_zLobatto);
            n = self.nEVP;
            self.xiLobatto = (Lxi/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + min(xi_zLobatto);
            
            % Now we need z on the \xi grid
            self.z_xiLobatto = interp1(xi_zLobatto, self.zLobatto, self.xiLobatto, 'spline');
            
            % and z_out on the \xi grid
            self.xiOut = interp1(self.zLobatto, xi_zLobatto, self.z, 'spline');
        end
        

        function self = InitializeWithFunction(self, rho, zMin, zMax, zOut)
            % Superclass calls this method upon initialization when it
            % determines that the input is given in functional form.
            %
            % The superclass will initialize zLobatto and rho_lobatto; this
            % class must initialize the sLobatto, z_sLobatto and sOut.
            %
            % Superclass initializes zLobatto and rho_lobatto
            InitializeWithFunction@InternalModesSpectral(self, rho, zMin, zMax, zOut);
            
            % Create the stretched grid \xi
            N_zLobatto = sqrt(InternalModesSpectral.ifct(self.N2_zCheb));
            self.N_zCheb = InternalModesSpectral.fct(N_zLobatto);
            xi_zLobatto = cumtrapz(self.zLobatto,N_zLobatto);
            Lxi = max(xi_zLobatto)-min(xi_zLobatto);
            n = self.nEVP;
            self.xiLobatto = (Lxi/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + min(xi_zLobatto);
            
            % Now we need z on the \xi grid
            self.z_xiLobatto = interp1(xi_zLobatto, self.zLobatto, self.xiLobatto, 'spline');
            
            % and z_out on the \xi grid
            self.xiOut = interp1(self.zLobatto, xi_zLobatto, self.z, 'spline');
        end
        
        function self = SetupEigenvalueProblem(self)
            % We will use the stretched grid to solve the eigenvalue
            % problem.
            self.xLobatto = self.xiLobatto;
            
            % The eigenvalue problem will be solved using N2 and N2z, so
            % now we need transformations to project them onto the
            % stretched grid
            T_zCheb_xiLobatto = InternalModesSpectral.ChebyshevTransformForGrid(self.zLobatto, self.z_xiLobatto);
            self.N2_xLobatto = T_zCheb_xiLobatto(self.N2_zCheb);
            self.Nz_xLobatto = T_zCheb_xiLobatto(self.Diff1_zCheb * self.N_zCheb);
            
            Lxi = max(self.xiLobatto) - min(self.xiLobatto);
            self.Diff1_xCheb = (2/Lxi)*InternalModesSpectral.ChebyshevDifferentiationMatrix( length(self.xiLobatto) );
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
            FOutFromGCheb = @(G_cheb,h) h * sqrt(self.N2) .* self.T_xCheb_zOut(self.Diff1_xCheb*G_cheb);
            GFromGCheb = @(G_cheb,h) InternalModesSpectral.ifct(G_cheb);
            FFromGCheb = @(G_cheb,h) h * sqrt(self.N2_xLobatto) .* InternalModesSpectral.ifct( self.Diff1_xCheb*G_cheb );
            GNorm = @(Gj) abs(sum(self.Int_xCheb .*InternalModesSpectral.fct((1/self.g) * (self.N2_xLobatto - self.f0*self.f0) .* ( self.N2_xLobatto.^(-0.5) ) .* Gj .^ 2)));
            FNorm = @(Fj) abs(sum(self.Int_xCheb .*InternalModesSpectral.fct((1/self.Lz) * (Fj.^ 2) .* ( self.N2_xLobatto.^(-0.5) ))));
            [F,G,h] = ModesFromGEP(self,A,B,hFromLambda,GFromGCheb,FFromGCheb,GNorm,FNorm,GOutFromGCheb,FOutFromGCheb);
        end
    end
    
end
