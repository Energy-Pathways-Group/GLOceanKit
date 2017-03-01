classdef InternalModesSpectral < InternalModesBase
    properties (Access = public)
        rho
        N2
    end
    
    properties (Dependent)
        rho_z
        rho_zz
    end
    
    % Need to establish notation:
    %
    % Dimensions ? like z ? have grids, which may differ. zIn is the
    % dimension z, on some grid given by the user. zLobatto is the
    % diension z, on a Lobatto grid. zCheb would be used for the Chebyshev
    % polynomial representation of a function defined on z. Note that the
    % base class property 'z' is really zOut.
    %
    % Functions ? if they're analytical, they're defined on a dimension,
    % not a grid. If they're gridded, then the grid must also be specified.
    % Thus, the function rho(z) when analytically defined will be given as
    % either just rho (when obvious) or rho_z when necessary. If it's
    % gridded, then expect rho_zIn or rho_zLobatto.
    %
    % Transformations ? a transformation often a matrix (or generally just
    % a function) that transformation a function from one basis to another.
    % A common scenario here is that we need to go from zLobatto to zCheb.
    % In our cases we're usually going from one gridded dimension to
    % another. Inverse Chebyshev transformations are denote by T, or more
    % specifically T_zCheb_zLobatto. Because T denotes Chebyshev, we will
    % simply reduce this to T_zLobatto
    %
    % The underscore in the case of rho_z is certainly recognized as a
    % derivative (because of LaTex notation), so it is awkward with the
    % notion used here, where the underscore is used to separate the
    % function/transformation name from the grid. Thus, we leave rho_z and
    % rho_zz as the only two exceptions to the aforementioned notation
    % because they are the only public facing interface.
    properties %(Access = private)
        zLobatto            % z_in coordinated, on Chebyshev extrema/Lobatto grid
        rho_zLobatto        % rho on the above zLobatto
        rho_zCheb           % rho in spectral space
        N2_zCheb            % N2 in spectral space
        N2_zLobatto         % N2 on the above zLobatto
        Diff1_zCheb         % single derivative in spectral space
        T_zLobatto, Tz_zLobatto, Tzz_zLobatto        % Chebyshev polys (and derivs) on the zLobatto
        
        doesOutputGridSpanDomain = 0    % if not, the output grid can't be used for normalization
        T_zCheb_zOut                          % *function* handle that transforms from spectral to z_out
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = InternalModesSpectral(rho, z_in, z_out, latitude)
            self@InternalModesBase(rho,z_in,z_out,latitude);
            
            self.rho_zCheb = fct(self.rho_zLobatto);
            self.Diff1_zCheb = (2/self.Lz)*ChebyshevDifferentiationMatrix( length(self.zLobatto) );
            self.N2_zCheb= -(self.g/self.rho0)*self.Diff1_zCheb*self.rho_zCheb;
            self.N2_zLobatto = ifct(self.N2_zCheb);
            
            [self.T_zLobatto,self.Tz_zLobatto,self.Tzz_zLobatto] = ChebyshevPolynomialsOnGrid( self.zLobatto, length(self.zLobatto) );
            [self.T_zCheb_zOut, self.doesOutputGridSpanDomain] = ChebyshevTransformForGrid(self.zLobatto, self.z);           
            
            self.rho = self.T_zCheb_zOut(self.rho_zCheb);
            self.N2 = self.T_zCheb_zOut(self.N2_zCheb);
        end
        
        % Superclass calls this method upon initialization when it
        % determines that the input is given in gridded form. The goal is
        % to initialize zLobatto and rho_zLobatto.
        function self = InitializeWithGrid(self, rho, z_in)
            if (z_in(2) - z_in(1)) > 0
                z_in = flip(z_in);
                rho = flip(rho);
            end
            
            if IsChebyshevGrid(z_in) == 1
                self.zLobatto = z_in;
                self.rho_zLobatto = rho;
            else
                self.zLobatto = FindSmallestChebyshevGridWithNoGaps( z_in ); % z, on a chebyshev grid
                self.rho_zLobatto = interp1(z_in, rho, self.zLobatto, 'spline'); % rho, interpolated to that grid
            end
        end
        
        % Superclass calls this method upon initialization when it
        % determines that the input is given in functional form. The goal
        % is to initialize zLobatto and rho_zLobatto.
        function self = InitializeWithFunction(self, rho, z_min, z_max, z_out)
            % We have an analytical function describing density, so lets
            % place points on a Chebyshev grid such that z_out is fully
            % resolved.
            if (z_out(2) - z_out(1)) > 0 % make z_out decreasing
                z_out = flip(z_out);
            end
            
            z_in = z_out;
            % Note that z_out might not span the whole domain, so we may
            % need to add z_min and z_max to the end points.
            if z_max > z_in(1)
                z_in = cat(1,z_max,z_in);
            end
            if z_min < z_in(end)
                z_in = cat(1,z_in,z_min);
            end
            
            if IsChebyshevGrid(z_in) == 1
                self.zLobatto = z_in;
            else
                self.zLobatto = FindSmallestChebyshevGridWithNoGaps( z_in ); % z, on a chebyshev grid
            end
            
            self.rho_zLobatto = rho(self.zLobatto);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computed (dependent) properties
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function value = get.rho_z(self)
            value = self.T_zCheb_zOut(self.Diff1_zCheb * self.rho_zCheb);
        end
        
        function value = get.rho_zz(self)
            value = self.T_zCheb_zOut(self.Diff1_zCheb * self.Diff1_zCheb * self.rho_zCheb);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computation of the modes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [F,G,h] = ModesAtWavenumber(self, k )
            % The eigenvalue equation is,
            % G_{zz} - K^2 G = \frac{f_0^2 -N^2}{gh_j}G
            % A = \frac{g}{f_0^2 -N^2} \left( \partial_{zz} - K^2*I \right)
            % B = I
            T = self.T_zLobatto;
            Tz = self.Tz_zLobatto;
            Tzz = self.Tzz_zLobatto;
            n = length(self.zLobatto);
            
            A = (Tzz - k*k*eye(n)*T);
            B = diag((self.f0*self.f0-self.N2_zLobatto)/self.g)*T;
            
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
            GFromGCheb = @(G_cheb,h) self.T_zCheb_zOut(G_cheb);
            FFromGCheb = @(G_cheb,h) h * self.T_zCheb_zOut(self.Diff1_zCheb*G_cheb);
            [F,G,h] = ModesFromGEP(self,A,B,hFromLambda,GFromGCheb,FFromGCheb);
        end
        
        % Take matrices A and B from the generalized eigenvalue problem
        % (GEP) and returns F,G,h. The h_func parameter is a function that
        % returns the eigendepth, h, given eigenvalue lambda from the GEP.
        function [F,G,h] = ModesFromGEP(self,A,B,hFromLambda,GFromGCheb,FFromGCheb)
            [V,D] = eig( A, B );
            
            [lambda, permutation] = sort(abs(diag(D)),1,'ascend');
            G_cheb=V(:,permutation);
            h = hFromLambda(lambda.');
            
            n = size(G_cheb,1);
            F = zeros(length(self.z),n);
            G = zeros(length(self.z),n);
            
            if self.doesOutputGridSpanDomain == 1
                for j=1:n
                    G(:,j) = GFromGCheb(G_cheb(:,j),h(j));
                    F(:,j) = FFromGCheb(G_cheb(:,j),h(j));
                end
                [F,G] = self.NormalizeModes(F,G,self.z);
            else
                error('This normalization condition is not yet implemented!')
            end
        end
        

        
    end
    
end


