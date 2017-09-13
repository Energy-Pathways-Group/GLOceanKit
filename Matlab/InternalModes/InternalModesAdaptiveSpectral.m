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
        N2_zLobatto         % Needs to be cached, because its used each time we create a new grid
        xi_zLobatto 
        
        gridFrequency        % The value of omega used for s = int sqrt(N^2 - omega^2) dz
        xiLobatto            % stretched density coordinate, on Chebyshev extrema/Lobatto grid
        z_xiLobatto          % The value of z, at the sLobatto points
        xiOut                % desired locations of the output in s-coordinate (deduced from z_out)
        
        N_zCheb
        Nz_xLobatto     	% (d/dz)N on the xiLobatto grid
        
        zBoundaries                 % z-location of the boundaries (end points plus turning points).
        xiBoundaries                % xi-location of the boundaries (end points plus turning points).
        nEquations
        boundaryIndicesStart        % indices of the boundaries into xiLobatto
        boundaryIndicesEnd          % indices of the boundaries into xiLobatto
        Lxi                         % array of length(nEquations) with the length of each EVP domain in xi coordinates
        
        T_xCheb_zOut_Transforms     % cell array containing function handles
        T_xCheb_zOut_fromIndices    % cell array with indices into the xLobatto grid
        T_xCheb_zOut_toIndices      % cell array with indices into the xiOut grid
        

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
            self.CreateGridForFrequency(omega);
            
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
            
            % now couple the equations together
            for i=2:self.nEquations
                eq1indices = self.boundaryIndicesStart(i-1):self.boundaryIndicesEnd(i-1);
                eq2indices = self.boundaryIndicesStart(i):self.boundaryIndicesEnd(i);
                         
                n = self.boundaryIndicesEnd(i-1);
                m = self.boundaryIndicesStart(i);
                % continuity in f
                
                A(n,eq1indices) = T(n,eq1indices);
                A(n,eq2indices) = -T(m,eq2indices);
                B(n,:) = 0;
                
                % continuity in df/dx
                A(m,eq1indices) = Tz(n,eq1indices);
                A(m,eq2indices) = -Tz(m,eq2indices);
                B(m,:) = 0;
            end
            
            [F,G,h] = self.ModesFromGEPWKBSpectral(A,B);
        end
 
        function v_xCheb = T_xLobatto_xCheb( self, v_xLobatto)
            % transform from xLobatto basis to xCheb basis
            v_xCheb = zeros(size(self.xiLobatto));
            for i=1:self.nEquations
                indices = self.boundaryIndicesStart(i):self.boundaryIndicesEnd(i);
                v_xCheb(indices) = InternalModesSpectral.fct(v_xLobatto(indices));
            end
        end
        
        function v_xLobatto = T_xCheb_xLobatto( self, v_xCheb)
            % transform from xCheb basis to xLobatto
            v_xLobatto = zeros(size(self.xiLobatto));
            for i=1:self.nEquations
                indices = self.boundaryIndicesStart(i):self.boundaryIndicesEnd(i);
                v_xLobatto(indices) = InternalModesSpectral.ifct(v_xCheb(indices));
            end
        end
        
        function v_zOut = T_xCheb_zOutFunction( self, v_xCheb )
            % transform from xCheb basis to zOut
            v_zOut = zeros(size(self.xiOut));
            for i = 1:length(self.T_xCheb_zOut_Transforms)
                T = self.T_xCheb_zOut_Transforms{i};
                v_zOut(self.T_xCheb_zOut_toIndices{i}) = T( v_xCheb(self.T_xCheb_zOut_fromIndices{i}) );
            end
        end
        
        function vx = Diff1_xChebFunction( self, v )
            % differentiate a vector in the compound Chebyshev xi basis
            vx = zeros(size(v));
            for iEquation = 1:self.nEquations
                indices = self.boundaryIndicesStart(iEquation):self.boundaryIndicesEnd(iEquation);
                vx(indices) = (2/self.Lxi(iEquation))*InternalModesSpectral.DifferentiateChebyshevVector( v(indices) );
            end
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
            
            self.InitializeStretchedCoordinates()
            self.CreateGridForFrequency(0);
        end
        

        function self = InitializeWithFunction(self, rho, zMin, zMax, zOut)
            % Superclass calls this method upon initialization when it
            % determines that the input is given in functional form.
            %
            % The superclass will initialize zLobatto and rho_lobatto; this
            % class must initialize the sLobatto, z_sLobatto and sOut.
            InitializeWithFunction@InternalModesSpectral(self, rho, zMin, zMax, zOut);
            
            self.InitializeStretchedCoordinates()
            self.CreateGridForFrequency(self.f0);
        end
        
        function InitializeStretchedCoordinates(self)
            self.N2_zLobatto = InternalModesSpectral.ifct(self.N2_zCheb);
            N_zLobatto = sqrt(self.N2_zLobatto);
            self.N_zCheb = InternalModesSpectral.fct(N_zLobatto);
            self.xi_zLobatto = cumtrapz(self.zLobatto,N_zLobatto);
            
            % Now we need z on the \xi grid
            self.z_xiLobatto = interp1(self.xi_zLobatto, self.zLobatto, self.xiLobatto, 'spline');
            
            % and z_out on the \xi grid
            self.xiOut = interp1(self.zLobatto, self.xi_zLobatto, self.z, 'spline');
        end
        
        function CreateGridForFrequency(self,omega)
            self.gridFrequency = omega;
            
            % Now find (roughly) the turning points, if any
            N2Omega2_zLobatto = self.N2_zLobatto - self.gridFrequency*self.gridFrequency;
            a = N2Omega2_zLobatto; a(a>=0) = 1; a(a<0) = 0;
            turningIndices = find(diff(a)~=0);
            nTP = length(turningIndices);
            zTP = zeros(nTP,1);
            for i=1:nTP
                fun = @(z) interp1(self.zLobatto,N2Omega2_zLobatto,z,'spline');
                zTP(i) = fzero(fun,self.zLobatto(turningIndices(i)));
            end  
            self.zBoundaries = [self.zLobatto(1); zTP; self.zLobatto(end)];
            self.xiBoundaries = interp1(self.zLobatto,self.xi_zLobatto,self.zBoundaries,'spline');
            
            % We will be coupling nTP+1 EVPs together. We need to
            % distribute the user requested points to each of these EVPs.
            % For this first draft, we simply evenly distribute the points.
            self.nEquations = nTP+1;
            nPoints = floor(self.nEVP/self.nEquations);
            nEVPPoints = nPoints*ones(self.nEquations,1);
            nEVPPoints(end) = nEVPPoints(end) + self.nEVP - nPoints*self.nEquations; % add any extra points to the end
            % A boundary point is repeated at the start of each EVP
            self.boundaryIndicesStart = cumsum( [1; nEVPPoints(1:end-1)] );
            self.boundaryIndicesEnd = self.boundaryIndicesStart + nEVPPoints-1;
            
            % Now we walk through the equations, and create a lobatto grid
            % for each equation.
            self.xiLobatto = zeros(self.nEVP,1);
            self.Int_xCheb = zeros(self.nEVP,1);
            self.T_xLobatto = zeros(self.nEVP,self.nEVP);
            self.Tx_xLobatto = zeros(self.nEVP,self.nEVP);
            self.Txx_xLobatto = zeros(self.nEVP,self.nEVP);
            for i=1:self.nEquations
                self.Lxi(i) = max(self.xiBoundaries(i+1),self.xiBoundaries(i)) - min(self.xiBoundaries(i+1),self.xiBoundaries(i));
                n = nEVPPoints(i);
                indices = self.boundaryIndicesStart(i):self.boundaryIndicesEnd(i);
                xLobatto = (self.Lxi(i)/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + min(self.xiBoundaries(i+1),self.xiBoundaries(i));
                self.xiLobatto(indices) = xLobatto;
                
                [T,Tx,Txx] = InternalModesSpectral.ChebyshevPolynomialsOnGrid( xLobatto, length(xLobatto) );
                self.T_xLobatto(indices,indices) = T;
                self.Tx_xLobatto(indices,indices) = Tx;
                self.Txx_xLobatto(indices,indices) = Txx;
                
                % We use that \int_{-1}^1 T_n(x) dx = \frac{(-1)^n + 1}{1-n^2}
                % for all n, except n=1, where the integral is zero.
                np = (0:(n-1))';
                Int = -(1+(-1).^np)./(np.*np-1);
                Int(2) = 0;
                Int = self.Lxi(i)/2*Int;
                self.Int_xCheb(indices) = Int;
            end
            
            % Now we need z on the \xi grid
            self.z_xiLobatto = interp1(self.xi_zLobatto, self.zLobatto, self.xiLobatto, 'spline');
            
            % and z_out on the \xi grid
            self.xiOut = interp1(self.zLobatto, self.xi_zLobatto, self.z, 'spline');
            
            % We will use the stretched grid to solve the eigenvalue
            % problem.
            self.xLobatto = self.xiLobatto;
            
            % The eigenvalue problem will be solved using N2 and N2z, so
            % now we need transformations to project them onto the
            % stretched grid
            T_zCheb_xiLobatto = InternalModesSpectral.ChebyshevTransformForGrid(self.zLobatto, self.z_xiLobatto);
            self.N2_xLobatto = T_zCheb_xiLobatto(self.N2_zCheb);
            self.Nz_xLobatto = T_zCheb_xiLobatto(self.Diff1_zCheb(self.N_zCheb));
            
            self.T_xCheb_zOut_Transforms = cell(0,0);
            self.T_xCheb_zOut_fromIndices = cell(0,0);
            self.T_xCheb_zOut_toIndices = cell(0,0);
            for i=1:self.nEquations
                if i == self.nEquations % upper and lower boundary included
                    toIndices = find( self.xiOut <= self.xiBoundaries(i) & self.xiOut >= self.xiBoundaries(i+1) );
                else % upper boundary included, lower boundary excluded (default)
                    toIndices = find( self.xiOut <= self.xiBoundaries(i) & self.xiOut > self.xiBoundaries(i+1) );
                end
                if ~isempty(toIndices)
                    index = length(self.T_xCheb_zOut_Transforms)+1;
                    fromIndices = self.boundaryIndicesStart(i):self.boundaryIndicesEnd(i);
                    
                    self.T_xCheb_zOut_fromIndices{index} = fromIndices;
                    self.T_xCheb_zOut_toIndices{index} = toIndices;
                    self.T_xCheb_zOut_Transforms{index} = InternalModesSpectral.ChebyshevTransformForGrid(self.xiLobatto(fromIndices), self.xiOut(toIndices));
                end
            end
            
            self.T_xCheb_zOut = @(v) self.T_xCheb_zOutFunction(v);
            self.Diff1_xCheb = @(v) self.Diff1_xChebFunction(v);
        end
                      
        function self = SetupEigenvalueProblem(self)            

            
        end
    end
    
    methods (Access = private)             
        function [F,G,h] = ModesFromGEPWKBSpectral(self,A,B)
            % This function is an intermediary used by ModesAtFrequency and
            % ModesAtWavenumber to establish the various norm functions.
            hFromLambda = @(lambda) 1.0 ./ lambda;
            GOutFromGCheb = @(G_cheb,h) self.T_xCheb_zOut(G_cheb);
            FOutFromGCheb = @(G_cheb,h) h * sqrt(self.N2) .* self.T_xCheb_zOut(self.Diff1_xChebFunction(G_cheb));
            GFromGCheb = @(G_cheb,h) self.T_xCheb_xLobatto(G_cheb);
            FFromGCheb = @(G_cheb,h) h * sqrt(self.N2_xLobatto) .* self.T_xCheb_xLobatto(self.Diff1_xChebFunction(G_cheb));
            GNorm = @(Gj) abs(sum(self.Int_xCheb .* self.T_xLobatto_xCheb((1/self.g) * (self.N2_xLobatto - self.f0*self.f0) .* ( self.N2_xLobatto.^(-0.5) ) .* Gj .^ 2)));
            FNorm = @(Fj) abs(sum(self.Int_xCheb .* self.T_xLobatto_xCheb((1/self.Lz) * (Fj.^ 2) .* ( self.N2_xLobatto.^(-0.5) ))));
            [F,G,h] = ModesFromGEP(self,A,B,hFromLambda,GFromGCheb,FFromGCheb,GNorm,FNorm,GOutFromGCheb,FOutFromGCheb);
        end
    end
    
end
