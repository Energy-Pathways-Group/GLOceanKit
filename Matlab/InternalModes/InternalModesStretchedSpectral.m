classdef InternalModesStretchedSpectral < InternalModesSpectral
    properties %(Access = private)    
        sLobatto            % stretched density coordinate, on Chebyshev extrema/Lobatto grid
        z_sLobatto          % The value of z, at the sLobatto points
        T_zCheb_sLobatto    % Matrix that transforms from zCheb to z_sLobatto (or sLobatto)
        T_sLobatto, Ts_sLobatto, Tss_sLobatto        % Chebyshev polys (and derivs) on the zLobatto
        Diff1_sCheb         % single derivative in spectral space
        sOut               % desired locations of the output in s-coordinate (deduced from z_out)
        
        N2_sLobatto     	% N2 on the z_sLobatto grid
        N2z_sLobatto    	% (d/dz)N2 on the z_sLobatto grid
        
        T_sCheb_zOut        %
    end
    

    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = InternalModesStretchedSpectral(rho, z_in, z_out, latitude, varargin)
            self@InternalModesSpectral(rho,z_in,z_out,latitude, varargin{:});
            
            Ls = max(self.sLobatto)-min(self.sLobatto);
            self.Diff1_sCheb = (2/Ls)*ChebyshevDifferentiationMatrix( length(self.sLobatto) );
            [self.T_sLobatto,self.Ts_sLobatto,self.Tss_sLobatto] = ChebyshevPolynomialsOnGrid( self.sLobatto, length(self.sLobatto) );
            [self.T_sCheb_zOut, self.doesOutputGridSpanDomain] = ChebyshevTransformForGrid(self.sLobatto, self.sOut);
            
            
            % The eigenvalue problem will be solved using N2 and N2z, so
            % now we need transformations to project them onto the
            % stretched grid
            t = acos((2/self.Lz)*(self.z_sLobatto-min(self.zLobatto)) - 1);
            self.T_zCheb_sLobatto = zeros(length(t),length(self.zLobatto));
            for iPoly=0:(length(self.zLobatto)-1)
                self.T_zCheb_sLobatto(:,iPoly+1) = cos(iPoly*t);
            end
            
            self.N2_sLobatto = self.T_zCheb_sLobatto * (self.N2_zCheb);
            self.N2z_sLobatto = self.T_zCheb_sLobatto * (self.Diff1_zCheb * self.N2_zCheb);
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
            
            s = -self.g*rho/self.rho0;
            L_s = max(s)-min(s);
%             warning('This choice of N_points is not good.')
            N_points = 2*length(z_in); % More is better, but this should be sufficient.
            xi=(0:N_points-1)';
            self.sLobatto = abs(L_s/2)*(cos(xi*pi/(N_points-1))+1) + min(s); % s Chebyshev grid on the s-coordinate
            self.z_sLobatto = interp1(s, z_in, self.sLobatto, 'spline'); % z, evaluated on that s grid
            
            self.sOut = interp1(z_in, s, self.z, 'spline'); % z, evaluated on that s grid
        end
        
        % Superclass calls this method upon initialization when it
        % determines that the input is given in functional form.
        %
        % The superclass will initialize zLobatto and rho_lobatto;
        % this class must initialize the sLobatto, z_sLobatto and
        % sOut.
        function self = InitializeWithFunction(self, rho, z_min, z_max, z_out)
            % Superclass initializes zLobatto and rho_lobatto
            InitializeWithFunction@InternalModesSpectral(self, rho, z_min, z_max, z_out);
            
            % Create a stretched grid that includes all the points of z_out
            s = @(z) -self.g*rho(z)/self.rho0;
            s_min = s(z_max);
            s_max = s(z_min);
            
            s_grid = s(z_out);
            if (s_grid(2) - s_grid(1)) > 0 % make s_grid decreasing
                s_grid = flip(s_grid);
            end
            if s_max > s_grid(1)
                s_grid = cat(1,s_max,s_grid);
            end
            if s_min < s_grid(end)
                s_grid = cat(1,s_grid,s_min);
            end
            
            if IsChebyshevGrid(s_grid) == 1
                self.sLobatto = s_grid;
            else
                % This will likely be waaay overkill in most cases... need
                % smarter logic here.
                self.sLobatto = FindSmallestChebyshevGridWithNoGaps(s_grid); % z, on a chebyshev grid
            end
            
            % Now create a transformation for functions defined on
            % z_lobatto to (spectrally) take them into s_lobatto.
            % We use the fact that we have a function handle to iteratively
            % improve this projection.
            self.z_sLobatto = interp1(s(self.zLobatto), self.zLobatto, self.sLobatto, 'spline');
%             for i=1:5
%                 self.z_sLobatto = interp1(s(self.z_sLobatto), self.z_sLobatto, self.sLobatto, 'spline');
%                 fprintf('max diff: %f\n', max(s(self.z_sLobatto) - self.sLobatto)/max(self.sLobatto));
%             end
            
            self.sOut = s(z_out);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computation of the modes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [F,G,h] = ModesAtWavenumber(self, k )
            T = self.T_sLobatto;
            Tz = self.Ts_sLobatto;
            Tzz = self.Tss_sLobatto;
            n = length(self.sLobatto);
            
            A = diag(self.N2_sLobatto .* self.N2_sLobatto)*Tzz + diag(self.N2z_sLobatto)*Tz - k*k*T;
            B = diag( (self.f0*self.f0 - self.N2_sLobatto)/self.g )*T;
            
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
            GFromGCheb = @(G_cheb,h) self.T_sCheb_zOut(G_cheb);
            FFromGCheb = @(G_cheb,h) h * self.N2 .* self.T_sCheb_zOut(self.Diff1_sCheb*G_cheb);
            [F,G,h] = self.ModesFromGEP(A,B,hFromLambda,GFromGCheb,FFromGCheb);
        end
        
        function [F,G,h] = ModesAtFrequency(self, omega )
            error('This function is not yet implemented!');
        end
        
    end
    
end
