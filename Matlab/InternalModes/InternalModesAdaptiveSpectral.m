classdef InternalModesAdaptiveSpectral < InternalModesWKBSpectral
    % This class solves the vertical eigenvalue problem on a WKB stretched
    % density coordinate grid using Chebyshev polynomials.
    %
    % See InternalModesBase for basic usage information.
    %
    % This class uses the coordinate
    %   xi = \int_{-Lz}^0 \sqrt(-(g/rho0)*rho_z) dz
    % to solve the EVP.
    %
    % This uses the same WKB coordinates as its superclass. The difference
    % is that the collocation points within that coordinate may be
    % different when more than one equation is used. Specifically, the
    % xLobatto grid is not actually Gauss-Lobatto.
    %
    %   See also INTERNALMODES, INTERNALMODESBASE, INTERNALMODESSPECTRAL,
    %   INTERNALMODESDENSITYSPECTRAL, and INTERNALMODESFINITEDIFFERENCE.
    %
    %   Jeffrey J. Early
    %   jeffrey@jeffreyearly.com
    %
    %   March 14th, 2017        Version 1.0
    
    properties %(Access = private)
        x_zLobatto                  % x (xi) coordinate on the zLobatto grid                
        N_zCheb
        
        zBoundaries                 % z-location of the boundaries (end points plus turning points).
        xiBoundaries                % xi-location of the boundaries (end points plus turning points).
        nEquations
        Lxi                         % array of length(nEquations) with the length of each EVP domain in xi coordinates
        
        eqIndices                   % cell array containing indices into the *rows* for a given EVP
        polyIndices                 % cell array containing indices into the *coumns* for a given EVP
        xiIndices
        
        T_xCheb_zOut_Transforms     % cell array containing function handles
        T_xCheb_zOut_fromIndices    % cell array with indices into the xLobatto grid
        T_xCheb_zOut_toIndices      % cell array with indices into the xOut grid
    end
    
    methods  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = InternalModesAdaptiveSpectral(rho, z_in, z_out, latitude, varargin)
            self@InternalModesWKBSpectral(rho,z_in,z_out,latitude, varargin{:});
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computation of the modes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [F,G,h,omega,varargout] = ModesAtWavenumber(self, k, varargin )
            % We just need to make sure we're using the right grid,
            % otherwise we can use the superclass method as is.
            self.CreateGridForFrequency(0.0);
            if isempty(varargin)
                [F,G,h,omega] = ModesAtWavenumber@InternalModesSpectral(self,k);
            else
                varargout = cell(size(varargin));
                [F,G,h,omega,varargout{:}] = ModesAtWavenumber@InternalModesSpectral(self,k, varargin{:});
            end
        end

        function [F,G,h,k,varargout] = ModesAtFrequency(self, omega, varargin )
            % We just need to make sure we're using the right grid,
            % otherwise we can use the superclass method as is.
            self.CreateGridForFrequency(omega);
            if isempty(varargin)
                [F,G,h,k] = ModesAtFrequency@InternalModesSpectral(self,omega);
            else
                varargout = cell(size(varargin));
                [F,G,h,k,varargout{:}] = ModesAtFrequency@InternalModesSpectral(self,omega, varargin{:});
            end
        end

        function [A,B] = EigenmatricesForFrequency(self, omega )
            T = self.T_xLobatto;
            Tz = self.Tx_xLobatto;
            Tzz = self.Txx_xLobatto;
            n = self.nEVP;

            N2 = zeros(self.nEVP,1);
            Nz = zeros(self.nEVP,1);
            for i=1:self.nEquations
                N2(self.eqIndices{i}) = self.N2_xLobatto(self.xiIndices{i});
                Nz(self.eqIndices{i}) = self.Nz_xLobatto(self.xiIndices{i});
            end
            
            A = diag(N2)*Tzz + diag(Nz)*Tz;
            B = diag( (omega*omega - N2)/self.g )*T;
            
            [A,B] = self.ApplyBoundaryConditions(A,B);
            
            % now couple the equations together, using the gaps we left in
            % the matrices.
            for i=2:self.nEquations                
                % continuity in f      
                A(max(self.eqIndices{i-1})+1,self.polyIndices{i-1}) = T(max(self.eqIndices{i-1}),self.polyIndices{i-1});
                A(max(self.eqIndices{i-1})+1,self.polyIndices{i}) = -T(min(self.eqIndices{i}),self.polyIndices{i});
                B(max(self.eqIndices{i-1})+1,:) = 0;
                
                % continuity in df/dx
                A(max(self.eqIndices{i-1})+2,self.polyIndices{i-1}) = Tz(max(self.eqIndices{i-1}),self.polyIndices{i-1});
                A(max(self.eqIndices{i-1})+2,self.polyIndices{i}) = -Tz(min(self.eqIndices{i}),self.polyIndices{i});
                B(max(self.eqIndices{i-1})+2,:) = 0;
            end
        end
 
        function v_xCheb = T_xLobatto_xCheb( self, v_xLobatto)
            % transform from xLobatto basis to xCheb basis
            v_xCheb = zeros(self.nEVP,1);
            for i=1:self.nEquations
                v_xCheb(self.polyIndices{i}(1:length(self.xiIndices{i}))) = InternalModesSpectral.fct(v_xLobatto(self.xiIndices{i}));
            end
        end
        
        function v_xLobatto = T_xCheb_xLobatto( self, v_xCheb)
            % transform from xCheb basis to xLobatto
            v_xLobatto = zeros(size(self.xLobatto));
            for i=1:self.nEquations
                v_xLobatto(self.xiIndices{i}) = InternalModesSpectral.ifct(v_xCheb(self.polyIndices{i}(1:length(self.xiIndices{i}))));
            end
        end
        
        function v_zOut = T_xCheb_zOutFunction( self, v_xCheb )
            % transform from xCheb basis to zOut
            v_zOut = zeros(size(self.xOut));
            for i = 1:length(self.T_xCheb_zOut_Transforms)
                T = self.T_xCheb_zOut_Transforms{i};
                v_zOut(self.T_xCheb_zOut_toIndices{i}) = T( v_xCheb(self.T_xCheb_zOut_fromIndices{i}) );
            end
        end
        
        function vx = Diff1_xChebFunction( self, v )
            % differentiate a vector in the compound Chebyshev xi basis
            vx = zeros(size(v));
            for iEquation = 1:self.nEquations
                vx(self.polyIndices{iEquation}) = (2/self.Lxi(iEquation))*InternalModesSpectral.DifferentiateChebyshevVector( v(self.polyIndices{iEquation}) );
            end
        end
        
    end
    
    methods (Access = protected)
        
        function self = SetupEigenvalueProblem(self)
            SetupEigenvalueProblem@InternalModesWKBSpectral(self);

            self.CreateGridForFrequency(0.0);
        end
                               
        function CreateGridForFrequency(self,omega)
            if self.gridFrequency == omega
                return
            end
            previousNEquations = self.nEquations;
            
            self.gridFrequency = omega;            
            
            requiredDecay = 1e-5; % 1e-6 performs worse at high frequencies, 1e-4 performs worse at low frequencies. 1e-5 wins.
            
            % Just to keep things simple, we do our calculations on gridded
            % data, rather than working with the chebfun and bspline
            % objects.
            n = 2^15 + 1;
            zLobatto = (self.Lz/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + min(self.zDomain);
            N2_zLobatto = self.N2_function(zLobatto);
            [zBoundariesAndTPs, thesign] = InternalModesAdaptiveSpectral.FindWKBSolutionBoundaries(N2_zLobatto, zLobatto, omega, requiredDecay);
%           [zBoundariesAndTPs, thesign] = InternalModesAdaptiveSpectral.FindWKBSolutionBoundaries(self.N2_zLobatto, self.zLobatto, omega, requiredDecay);
%           [zBoundariesAndTPs2, thesign2] = self.FindWKBSolutionBoundariesSpectral(omega,requiredDecay);
            boundaries = self.x_function(zBoundariesAndTPs);
            
            % Distribute the allowed points in these regions, but remove
            % regions that end up with too few points.
            nEVPPoints = InternalModesAdaptiveSpectral.DistributePointsInRegionsWithMinimum( self.nEVP, boundaries, thesign );
            [minPoints, minIndex] = min(nEVPPoints);
            while ( minPoints < 5 )
                [boundaries, thesign] = InternalModesAdaptiveSpectral.RemoveRegionAtIndex( boundaries, thesign, minIndex );
                nEVPPoints = InternalModesAdaptiveSpectral.DistributePointsInRegionsWithMinimum( self.nEVP, boundaries, thesign );
                [minPoints, minIndex] = min(nEVPPoints);
            end
                        
            self.SetupCoupledEquationsAtBoundaries( boundaries, nEVPPoints );
            
            if self.nEquations ~= previousNEquations
                if self.nEquations == 1
                    fprintf(' The eigenvalue problem will be solved with %d points.\n', length(self.xLobatto));
                else
                    fprintf(' Separating the EVP into %d coupled EVPs with ',self.nEquations);
                    for i=1:length(nEVPPoints)
                        if i == length(nEVPPoints)
                            fprintf('%d',nEVPPoints(i));
                        else
                            fprintf('%d+',nEVPPoints(i));
                        end
                    end
                    fprintf(' points.\n');
                end
            end

            self.hFromLambda = @(lambda) 1.0 ./ lambda;
            self.GOutFromVCheb = @(G_cheb,h) self.T_xCheb_zOut(G_cheb);
            self.FOutFromVCheb = @(G_cheb,h) h * sqrt(self.N2) .* self.T_xCheb_zOut(self.Diff1_xChebFunction(G_cheb));
            self.GFromVCheb = @(G_cheb,h) self.T_xCheb_xLobatto(G_cheb);
            self.FFromVCheb = @(G_cheb,h) h * sqrt(self.N2_xLobatto) .* self.T_xCheb_xLobatto(self.Diff1_xChebFunction(G_cheb));
            self.GNorm = @(Gj) abs(Gj(1)*Gj(1) + sum(self.Int_xCheb .* self.T_xLobatto_xCheb((1/self.g) * (self.N2_xLobatto - self.f0*self.f0) .* ( self.N2_xLobatto.^(-0.5) ) .* Gj .^ 2)));
            self.FNorm = @(Fj) abs(sum(self.Int_xCheb .* self.T_xLobatto_xCheb((1/self.Lz) * (Fj.^ 2) .* ( self.N2_xLobatto.^(-0.5) ))));
        end
        
        function SetupCoupledEquationsAtBoundaries( self, boundaries, nEVPPoints )
            self.xiBoundaries = boundaries;
            self.zBoundaries = InternalModesSpectral.fInverseBisection(self.x_function,self.xiBoundaries,self.zMin,self.zMax,1e-12);
            self.zBoundaries(self.zBoundaries>self.zMax) = self.zMax;
            self.zBoundaries(self.zBoundaries<self.zMin) = self.zMin;
            self.Lxi = abs(diff(boundaries));
            self.nEquations = length(self.xiBoundaries)-1;
            
            % A boundary point is repeated at the start of each EVP
            boundaryIndicesStart = cumsum( [1; nEVPPoints(1:end-1)] );
            boundaryIndicesEnd = boundaryIndicesStart + nEVPPoints-1;
            
            % The final matrices in the EVP will be square, with each
            % column representing a Chebyshev polynomial, and each row
            % representing an equation (usually the poly at that point).
            %
            % The polyIndices are the indices into the column for a given
            % equation. There are no gaps in the polyIndices, so that an
            % equation with n points will have n polynomials.
            %
            % The eqIndices are the indices into the rows for a given
            % equations. There are gaps here, as we leave an extra two
            % rows at each equation coupling so that we can enforcing
            % continuity of f an df/dx.
            %
            % The xiIndices map the unique points in xi into
            % self.xLobatto. Each equation index corresponds to a grid
            % point... this is how we avoid those gaps we created.
            self.eqIndices = cell(self.nEquations,1);
            self.polyIndices = cell(self.nEquations,1);
            self.xiIndices = cell(self.nEquations,1);
            endXiIndex = 0;
            for i=1:self.nEquations
                self.polyIndices{i} = boundaryIndicesStart(i):boundaryIndicesEnd(i);
                if i==1 && i==self.nEquations
                    self.eqIndices{i} = boundaryIndicesStart(i):boundaryIndicesEnd(i);
                elseif i==1
                    self.eqIndices{i} = boundaryIndicesStart(i):(boundaryIndicesEnd(i)-1);
                elseif i==self.nEquations
                    self.eqIndices{i} = (boundaryIndicesStart(i)+1):boundaryIndicesEnd(i);
                else
                    self.eqIndices{i} = (boundaryIndicesStart(i)+1):(boundaryIndicesEnd(i)-1);
                end
                startXiIndex = endXiIndex + 1;
                endXiIndex = startXiIndex + length(self.eqIndices{i})-1;
                self.xiIndices{i} = startXiIndex:endXiIndex;
            end
            nGridPoints = self.nEVP - 2*(self.nEquations-1);
                        
            % Now we walk through the equations, and create a lobatto grid
            % for each equation.
            self.xLobatto = zeros(nGridPoints,1);
            self.Int_xCheb = zeros(self.nEVP,1);
            self.T_xLobatto = zeros(self.nEVP,self.nEVP);
            self.Tx_xLobatto = zeros(self.nEVP,self.nEVP);
            self.Txx_xLobatto = zeros(self.nEVP,self.nEVP);
            for i=1:self.nEquations
                n = length(self.eqIndices{i});
                m = length(self.polyIndices{i});
                xLobatto_local = (self.Lxi(i)/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + min(self.xiBoundaries(i+1),self.xiBoundaries(i));
                self.xLobatto(self.xiIndices{i}) = xLobatto_local;
                
                [T,Tx,Txx] = InternalModesSpectral.ChebyshevPolynomialsOnGrid( xLobatto_local, m );
                self.T_xLobatto(self.eqIndices{i},self.polyIndices{i}) = T;
                self.Tx_xLobatto(self.eqIndices{i},self.polyIndices{i}) = Tx;
                self.Txx_xLobatto(self.eqIndices{i},self.polyIndices{i}) = Txx;
                
                % We use that \int_{-1}^1 T_n(x) dx = \frac{(-1)^n + 1}{1-n^2}
                % for all n, except n=1, where the integral is zero.
                np = (0:(m-1))';
                Int = -(1+(-1).^np)./(np.*np-1);
                Int(2) = 0;
                Int = self.Lxi(i)/2*Int;
                self.Int_xCheb(self.polyIndices{i}) = Int;
            end
            
            % Now we need z on the \xi grid
            self.z_xLobatto = InternalModesSpectral.fInverseBisection(self.x_function,self.xLobatto,self.zMin,self.zMax,1e-12);
            self.z_xLobatto(self.z_xLobatto>self.zMax) = self.zMax;
            self.z_xLobatto(self.z_xLobatto<self.zMin) = self.zMin;
                        
            % The eigenvalue problem will be solved using N2 and N2z, so
            % now we need transformations to project them onto the
            % stretched grid
            self.N2_xLobatto = self.N2_function(self.z_xLobatto);
            self.Nz_xLobatto = self.Nz_function(self.z_xLobatto);
            
            self.T_xCheb_zOut_Transforms = cell(0,0);
            self.T_xCheb_zOut_fromIndices = cell(0,0);
            self.T_xCheb_zOut_toIndices = cell(0,0);
            for i=1:self.nEquations
                if i == self.nEquations % upper and lower boundary included
                    toIndices = find( self.xOut <= self.xiBoundaries(i) & self.xOut >= self.xiBoundaries(i+1) );
                else % upper boundary included, lower boundary excluded (default)
                    toIndices = find( self.xOut <= self.xiBoundaries(i) & self.xOut > self.xiBoundaries(i+1) );
                end
                if ~isempty(toIndices)
                    index = length(self.T_xCheb_zOut_Transforms)+1;
                    
                    self.T_xCheb_zOut_fromIndices{index} = self.polyIndices{i}(1:length(self.xiIndices{i}));
                    self.T_xCheb_zOut_toIndices{index} = toIndices;
                    self.T_xCheb_zOut_Transforms{index} = InternalModesSpectral.ChebyshevTransformForGrid(self.xLobatto(self.xiIndices{i}), self.xOut(toIndices));
                end
            end
            
            self.T_xCheb_zOut = @(v) self.T_xCheb_zOutFunction(v);
            self.Diff1_xCheb = @(v) self.Diff1_xChebFunction(v);
        end
        
        function [zBoundariesAndTPs, thesign] = FindWKBSolutionBoundariesSpectral(self, omega, requiredDecay)
            N2Omega2_zCheb = self.N2_zCheb;
            N2Omega2_zCheb(1) = N2Omega2_zCheb(1) - omega*omega;
            
            % boundaries and turning points, ordered top to bottom
            zTPs = self.FindTurningPointsAtFrequency(omega);
            zBoundariesAndTPs = [max(self.zDomain); sort(zTPs,'descend'); min(self.zDomain)];
            
            % sign of the solution in those regions
            midZ = zBoundariesAndTPs(1:end-1) + diff(zBoundariesAndTPs)/2;
            thesign = sign(InternalModesSpectral.ValueOfFunctionAtPointOnGrid(midZ,self.zLobatto,N2Omega2_zCheb));

            % sum of the oscillatory (positive sign) regions
            L_osc = 0.0;
            for i = 1:length(thesign)
                if thesign(i) > 0
                    L_osc = L_osc + InternalModesSpectral.IntegrateChebyshevVectorWithLimits(N2Omega2_zCheb,self.zLobatto,zBoundariesAndTPs(i+1),zBoundariesAndTPs(i));
                end
            end
            
            q_zCheb = InternalModesSpectral.IntegrateChebyshevVector( InternalModesSpectral.fct(sqrt(abs(self.N2_zLobatto - omega*omega))) );
            q = InternalModesSpectral.ifct(q_zCheb);
            
        end
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Static Methods (that do not depend on self)
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Static)
        
%         function [zBoundariesAndTPs, thesign] = FindWKBSolutionBoundaries(N2, omega, requiredDecay)
%             % N2 is either a chebfun or a bspline. bspline was written to
%             % respond to the same functions as chebfun.
%             
%             % First just find the turning points
%             N2Omega2 = N2 - omega*omega;
%             zBoundariesAndTPs = [min(N2Omega2.domain); reshape(roots(N2Omega2),[],1); max(N2Omega2.domain)];
%             midZ = zBoundariesAndTPs(1:end-1) + diff(zBoundariesAndTPs)/2;
%             thesign = sign( N2Omega2(midZ) );
%             
%             % Now sum up the total oscillatory region
%             L_osc = 0.0;
%             for i = 1:length(thesign)
%                 if thesign(i) > 0
%                     indices = boundaryIndices(i+1):-1:boundaryIndices(i);
%                     L_osc = L_osc + trapz(z(indices), abs(N2Omega2(indices)).^(1/2));
%                 end
%             end
%         end
        
        function [zBoundariesAndTPs, thesign] = FindWKBSolutionBoundaries(N2, z, omega, requiredDecay)
            % Find the theoretical solution boundaries of the WKB solution.
            % The requiredDecay is <=1 and find the point at which the
            % solution has decayed below that value.
            [zBoundariesAndTPs, thesign, boundaryIndices] = InternalModesSpectral.FindTurningPointBoundariesAtFrequency(N2, z, omega);
            N2Omega2 = N2 - omega*omega;
            
            L_osc = 0.0;
            for i = 1:length(thesign)
                if thesign(i) > 0
                    indices = boundaryIndices(i+1):-1:boundaryIndices(i);
                    L_osc = L_osc + trapz(z(indices), abs(N2Omega2(indices)).^(1/2));
                end
            end
            
            indicesToRemove = [];
            q = cumtrapz(z, abs(N2Omega2).^(1/2));
            for iInteriorPoint = 2:(length(zBoundariesAndTPs)-1)
                % xi is just q, but where its now 0 at the turning point of interest.
                xi = abs(q-interp1(z,q,zBoundariesAndTPs(iInteriorPoint)));
                if thesign( iInteriorPoint-1) < 0
                    % negative (decay) to the left, positive (oscillatory) to the right
                    indices = boundaryIndices(iInteriorPoint):-1:boundaryIndices(iInteriorPoint-1);
                    decay = InternalModesAdaptiveSpectral.WKBDecaySolution(xi(indices), L_osc, N2Omega2(indices));
                    
                    decayIndex = find( decay/max(decay) < requiredDecay, 1, 'first');
                    if isempty(decayIndex) || indices(decayIndex) < boundaryIndices(iInteriorPoint-1)
                        indicesToRemove(end+1) = iInteriorPoint;
                    else
                        boundaryIndices(iInteriorPoint) = indices(decayIndex);
                    end
                else
                    % positive (oscillatory) to the left, negative (decay) to the right
                    indices = (boundaryIndices(iInteriorPoint)+1):boundaryIndices(iInteriorPoint+1);
                    decay = InternalModesAdaptiveSpectral.WKBDecaySolution(xi(indices), L_osc, N2Omega2(indices));
                    
                    decayIndex = find( decay/max(decay) < requiredDecay, 1, 'first');
                    if isempty(decayIndex) || indices(decayIndex) > boundaryIndices(iInteriorPoint+1)
                        indicesToRemove(end+1) = iInteriorPoint;
                    else
                        boundaryIndices(iInteriorPoint) = indices(decayIndex);
                    end
                end
                
            end
                    
            if ~isempty(indicesToRemove)
                zBoundariesAndTPs( indicesToRemove ) = [];
                boundaryIndices( indicesToRemove ) = [];
                % the sign will always alternate. -1/+1
                if min(indicesToRemove) == 2
                    thesign = -1*ones(1,length(boundaryIndices)-1);
                    thesign(1) = 1;
                    thesign = cumprod(thesign);
                else
                    thesign = thesign(1:(length(boundaryIndices)-1));
                    if length(thesign)>1
                        thesign(2:end) = -1;
                        thesign = cumprod(thesign);
                    end
                end
            end
            
            zBoundariesAndTPs = z(boundaryIndices); 
        end
        
        function decay = WKBDecaySolution(xi, L_osc, N2Omega2)
            % The decay part of the lowest F-mode. xi should be > 0
            c = L_osc./((3/4)*pi);
            q = xi * (1./c);
            eta = (3*abs(q)/2).^(2/3) ;
            eta_z = -(sqrt(abs(N2Omega2))*(1./c)) .* (3*abs(q)/2).^(-1/3);
            decay = (-eta./(N2Omega2)).^(1/4) .* eta_z .* airy(1,eta);
            
            % Here's the simpler, approximated solution that has problems
            % near 0.
            % decay = (abs(N2Omega2).^(1/4)).*exp( -(3/4)*pi*xi/L_osc);
        end
               
        function [newBoundaries, newSigns] = RemoveRegionAtIndex(oldBoundaries, oldSigns, index)
            % Given a list of boundaries (length n), and the signs of the
            % regions created by the boundaries (length n-1), this removes
            % the boundary at some index, and returns the appropriate signs
            newBoundaries = oldBoundaries;
            newSigns = oldSigns;
            if length(oldBoundaries) < 3
                return
            elseif index == 1
                newBoundaries(2) = [];
                newSigns(1) = [];
            elseif index == length(oldBoundaries)-1
                newBoundaries(length(oldBoundaries)-1) = [];
                newSigns(end) = [];
            else
                newBoundaries([index;index+1]) = [];
                newSigns([index-1;index]) = [];
            end
        end
                
        function nEVPPoints = DistributePointsInRegionsWithMinimum( nTotalPoints, boundaries, thesign )
            % This algorithm distributes the points/polynomials amongst the
            % different regions/equations
            L = abs(diff(boundaries));
            totalEquations = length(boundaries)-1;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This is a key parameter that can strongly influence the
            % results. We never let it drop below 6, and generally try to
            % keep it at 2^4 (1/16th) the total points).
            minPoints = max([round(2^(log2(nTotalPoints)-4))/sum(thesign<0) 6]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
            if totalEquations > 1
                nEVPPoints = zeros(totalEquations,1);
                
                nEVPPoints(thesign<0) = minPoints;
                
                nPositivePoints = nTotalPoints - minPoints*sum(thesign<0);
                indices = thesign>0;
                if ~isempty(indices)
                    relativeSize = sqrt( L(indices)./min(L(indices)));
                    nBase = nPositivePoints/sum( relativeSize );
                    nEVPPoints(indices) = floor( relativeSize * nBase );
                    nEVPPoints(indices(end)) = nEVPPoints(indices(end)) + nPositivePoints-sum(nEVPPoints(indices));
                end
                
            else
                nEVPPoints = nTotalPoints;
            end
        end
    end
    
end
