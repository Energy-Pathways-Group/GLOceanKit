%% InternalModesStretchedSpectral
% 
% This class solves the vertical eigenvalue problem on a stretched density
% coordinate grid using Chebyshev polynomials.
%
% Internally, sLobatto is the stretched density coordinate on a Chebyshev
% extrema/Lobatto grid. This is the grid upon which the eigenvalue problem
% is solved, and therefore the class uses the superclass properties denoted
% with 'x' instead of 's' when setting up the eigenvalue problem.
classdef InternalModesStretchedSpectral < InternalModesSpectral
    properties %(Access = private)    
        sLobatto            % stretched density coordinate, on Chebyshev extrema/Lobatto grid
        z_sLobatto          % The value of z, at the sLobatto points
        sOut                % desired locations of the output in s-coordinate (deduced from z_out)
        
        N2z_xLobatto    	% (d/dz)N2 on the z_sLobatto grid   
    end
    

    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = InternalModesStretchedSpectral(rho, z_in, z_out, latitude, varargin)
            self@InternalModesSpectral(rho,z_in,z_out,latitude, varargin{:});
        end
                
        % Superclass calls this method upon initialization when it
        % determines that the input is given in gridded form.
        %
        % The superclass will initialize zLobatto and rho_lobatto;
        % this class must initialize the sLobatto, z_sLobatto and
        % sOut.
        function self = InitializeWithGrid(self, rho, z_in)
            % Superclass initializes zLobatto and rho_lobatto
            InitializeWithGrid@InternalModesSpectral(self, rho, z_in);
            
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
            s = -self.g*rho/self.rho0;
            Ls = max(s)-min(s);
            self.sLobatto = (Ls/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + min(s);
            self.z_sLobatto = interp1(s, z_in, self.sLobatto, 'spline'); % z, evaluated on that s grid
            
            self.sOut = interp1(z_in, s, self.z, 'spline'); % z, evaluated on that s grid
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
            s = @(z) -self.g*rho(z)/self.rho0;
            
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
            sMin = s(zMax);
            sMax = s(zMin);
            Ls = sMax-sMin;
            self.sLobatto = (Ls/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + sMin;            
            
            % Now create a transformation for functions defined on
            % z_lobatto to (spectrally) take them into s_lobatto.
            % We use the fact that we have a function handle to iteratively
            % improve this projection.
            self.z_sLobatto = interp1(s(self.zLobatto), self.zLobatto, self.sLobatto, 'spline');
            for i=1:5
                self.z_sLobatto = interp1(s(self.z_sLobatto), self.z_sLobatto, self.sLobatto, 'spline');
                fprintf('max diff: %g\n', max(s(self.z_sLobatto) - self.sLobatto)/max(self.sLobatto));
            end
            
            self.sOut = s(zOut);
        end
        
        function self = SetupEigenvalueProblem(self)
            % We will use the stretched grid to solve the eigenvalue
            % problem.
            self.xLobatto = self.sLobatto;
            
            % The eigenvalue problem will be solved using N2 and N2z, so
            % now we need transformations to project them onto the
            % stretched grid
            T_zCheb_sLobatto = ChebyshevTransformForGrid(self.zLobatto, self.z_sLobatto);
            self.N2_xLobatto = T_zCheb_sLobatto(self.N2_zCheb);
            self.N2z_xLobatto = T_zCheb_sLobatto(self.Diff1_zCheb * self.N2_zCheb);
            
            Ls = max(self.sLobatto)-min(self.sLobatto);
            self.Diff1_xCheb = (2/Ls)*ChebyshevDifferentiationMatrix( length(self.sLobatto) );
            [self.T_xLobatto,self.Tx_xLobatto,self.Txx_xLobatto] = ChebyshevPolynomialsOnGrid( self.sLobatto, length(self.sLobatto) );
            [self.T_xCheb_zOut, self.doesOutputGridSpanDomain] = ChebyshevTransformForGrid(self.sLobatto, self.sOut);
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
            
%             a = diag(self.N2_xLobatto .* self.N2_xLobatto)*Tz;
%             b = diag(self.N2z_xLobatto)*Tz;
%             c = k*k*T;
%             d = diag( (self.f0*self.f0 - self.N2_xLobatto)/self.g )*T;
%             climits = [-10 -3];
%             figure
%             subplot(2,2,1)
%             pcolor(log10(abs(a))), shading flat, caxis(climits)
%             subplot(2,2,2)
%             pcolor(log10(abs(b))), shading flat, caxis(climits)
%             subplot(2,2,3)
%             pcolor(log10(abs(c))), shading flat, caxis(climits)
%             subplot(2,2,4)
%             pcolor(log10(abs(d))), shading flat, caxis(climits)
%             
            A = diag(self.N2_xLobatto .* self.N2_xLobatto)*Tzz + diag(self.N2z_xLobatto)*Tz - k*k*T;
            B = diag( (self.f0*self.f0 - self.N2_xLobatto)/self.g )*T;
            
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
            FFromGCheb = @(G_cheb,h) h * self.N2 .* self.T_xCheb_zOut(self.Diff1_xCheb*G_cheb);
            [F,G,h] = self.ModesFromGEP(A,B,hFromLambda,GFromGCheb,FFromGCheb);
        end
        
        function [F,G,h] = ModesAtFrequency(self, omega )
            error('This function is not yet implemented!');
        end
        
    end
    
end
