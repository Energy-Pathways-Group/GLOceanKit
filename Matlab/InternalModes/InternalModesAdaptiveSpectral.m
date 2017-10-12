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
        
        xiLobatto            % stretched density coordinate, on Chebyshev extrema/Lobatto grid
        z_xiLobatto          % The value of z, at the sLobatto points
        xiOut                % desired locations of the output in s-coordinate (deduced from z_out)
        
        N_zCheb
        Nz_xLobatto     	% (d/dz)N on the xiLobatto grid
        
        zBoundaries                 % z-location of the boundaries (end points plus turning points).
        xiBoundaries                % xi-location of the boundaries (end points plus turning points).
        nEquations
%         boundaryIndicesStart        % indices of the boundaries into xiLobatto
%         boundaryIndicesEnd          % indices of the boundaries into xiLobatto
        Lxi                         % array of length(nEquations) with the length of each EVP domain in xi coordinates
        
        eqIndices                   % cell array containing indices into the *rows* for a given EVP
        polyIndices                 % cell array containing indices into the *coumns* for a given EVP
        xiIndices
        
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
            self.CreateGridForFrequency(0.0);
            
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

            N2 = zeros(self.nEVP,1);
            Nz = zeros(self.nEVP,1);
            for i=1:self.nEquations
                N2(self.eqIndices{i}) = self.N2_xLobatto(self.xiIndices{i});
                Nz(self.eqIndices{i}) = self.Nz_xLobatto(self.xiIndices{i});
            end
            
            A = diag(N2)*Tzz + diag(Nz)*Tz;
            B = diag( (omega*omega - N2)/self.g )*T;
            
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
                % continuity in f      
                A(max(self.eqIndices{i-1})+1,self.polyIndices{i-1}) = T(max(self.eqIndices{i-1}),self.polyIndices{i-1});
                A(max(self.eqIndices{i-1})+1,self.polyIndices{i}) = -T(min(self.eqIndices{i}),self.polyIndices{i});
                B(max(self.eqIndices{i-1})+1,:) = 0;
                
                % continuity in df/dx
                A(max(self.eqIndices{i-1})+2,self.polyIndices{i-1}) = Tz(max(self.eqIndices{i-1}),self.polyIndices{i-1});
                A(max(self.eqIndices{i-1})+2,self.polyIndices{i}) = -Tz(min(self.eqIndices{i}),self.polyIndices{i});
                B(max(self.eqIndices{i-1})+2,:) = 0;
            end
            
            [F,G,h] = self.ModesFromGEPWKBSpectral(A,B);
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
            v_xLobatto = zeros(size(self.xiLobatto));
            for i=1:self.nEquations
                v_xLobatto(self.xiIndices{i}) = InternalModesSpectral.ifct(v_xCheb(self.polyIndices{i}(1:length(self.xiIndices{i}))));
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
                vx(self.polyIndices{iEquation}) = (2/self.Lxi(iEquation))*InternalModesSpectral.DifferentiateChebyshevVector( v(self.polyIndices{iEquation}) );
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
            self.CreateGridForFrequency(0.0);
        end
        

        function self = InitializeWithFunction(self, rho, zMin, zMax, zOut)
            % Superclass calls this method upon initialization when it
            % determines that the input is given in functional form.
            %
            % The superclass will initialize zLobatto and rho_lobatto; this
            % class must initialize the sLobatto, z_sLobatto and sOut.
            InitializeWithFunction@InternalModesSpectral(self, rho, zMin, zMax, zOut);
            
            self.InitializeStretchedCoordinates()
            self.CreateGridForFrequency(0.0);
        end
        
        function InitializeStretchedCoordinates(self)
            N_zLobatto = sqrt(self.N2_zLobatto);
            self.N_zCheb = InternalModesSpectral.fct(N_zLobatto);
            self.xi_zLobatto = cumtrapz(self.zLobatto,N_zLobatto);
            
            % Now we need z on the \xi grid
            self.z_xiLobatto = interp1(self.xi_zLobatto, self.zLobatto, self.xiLobatto, 'spline');
            
            % and z_out on the \xi grid
            self.xiOut = interp1(self.zLobatto, self.xi_zLobatto, self.z, 'spline');
        end
               
        function [zBoundariesAndTPs, thesign] = FindWKBSolutionBoundaries(self, omega)
            requiredDecay = 1e-6; % 1e-6 seems to prevent exp strat from doing any worse... going to 1e-5 will do slightly worse for some frequencies.
            
            [zBoundariesAndTPs, thesign, boundaryIndices] = self.FindTurningPointBoundariesAtFrequency(omega);
            
            L_osc = 0.0;
            for i = 1:length(thesign)
                if thesign(i) > 0
                    indices = boundaryIndices(i+1):-1:boundaryIndices(i);
                    L_osc = L_osc + trapz(self.zLobatto(indices), abs(N2Omega2_zLobatto(indices)).^(1/2));
                end
            end
            
            indicesToRemove = [];
            for iInteriorPoint = 2:(length(zBoundariesAndTPs)-1)
               if thesign( iInteriorPoint-1) < 0
                  % negative (decay) to the left, positive (oscillatory) to
                  % the right
                  indices = boundaryIndices(iInteriorPoint):-1:boundaryIndices(iInteriorPoint-1);
                  xi = cumtrapz(self.zLobatto(indices), abs(N2Omega2_zLobatto(indices)).^(1/2));
                  decay = (abs(N2Omega2_zLobatto(indices)).^(1/4)).*exp( -pi*xi/L_osc);
                  decayIndex = find( decay/max(decay) < requiredDecay, 1, 'first');
                  if isempty(decayIndex) || indices(decayIndex) < boundaryIndices(iInteriorPoint-1)
                      % need to remove this point
                      indicesToRemove(end+1) = iInteriorPoint;
                  else
                      boundaryIndices(iInteriorPoint) = indices(decayIndex);
                  end
               else
                   % positive (oscillatory) to the left, negative (decay) to
                   % the right
                   indices = (boundaryIndices(iInteriorPoint)+1):boundaryIndices(iInteriorPoint+1);
                   xi = -cumtrapz(self.zLobatto(indices), abs(N2Omega2_zLobatto(indices)).^(1/2));
                   decay = (abs(N2Omega2_zLobatto(indices)).^(1/4)).*exp( -pi*xi/L_osc);
                   decayIndex = find( decay/max(decay) < requiredDecay, 1, 'first');
                   if isempty(decayIndex) || indices(decayIndex) > boundaryIndices(iInteriorPoint+1)
                       % need to remove this point
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
            
            
            zBoundariesAndTPs = self.zLobatto(boundaryIndices);
            
        end
        
        function CreateGridForFrequency(self,omega)
            if self.gridFrequency == omega
                return
            end
            
            self.gridFrequency = omega;
                        
%             [zBoundariesAndTPs, thesign] = self.FindTurningPointBoundariesAtFrequency(self.gridFrequency);
%             xiBoundariesAndTPs = interp1(self.zLobatto,self.xi_zLobatto,zBoundariesAndTPs,'spline');
%             
%             % Extended the oscillatory regions by some predefined amount
%             % (L_osc)
%             [boundaries, thesign] = self.GrowOscillatoryRegions( xiBoundariesAndTPs, thesign );
            
            
            [zBoundariesAndTPs, thesign] = FindWKBSolutionBoundaries(self, omega);
            boundaries = interp1(self.zLobatto,self.xi_zLobatto,zBoundariesAndTPs,'spline');
            
            % Distribute the allowed points in these regions, but remove
            % regions that end up with too few points.
            nEVPPoints = self.DistributePointsInRegionsWithMinimum( self.nEVP, boundaries, thesign );
            [minPoints, minIndex] = min(nEVPPoints);
            while ( minPoints < 5 )
                [boundaries, thesign] = self.RemoveRegionAtIndex( boundaries, thesign, minIndex );
                nEVPPoints = self.DistributePointsInRegionsWithMinimum( self.nEVP, boundaries, thesign );
                [minPoints, minIndex] = min(nEVPPoints);
            end
                        
            self.SetupCoupledEquationsAtBoundaries( boundaries, nEVPPoints );
            
            nEVPPoints
%             self.Lxi
        end
        
        function [newBoundaries, newsigns] = GrowOscillatoryRegions(self, xiBoundariesAndTPs, thesign)
            LxiRegion = abs(diff(xiBoundariesAndTPs));
            newBoundaries = zeros(size(xiBoundariesAndTPs));
            newsigns = zeros(size(thesign));
            
            % We pick a length scale (in WKB coordinates) over which the
            % mode decays approximately three orders of magnitude.
            L_osc = 6*sum(LxiRegion(thesign>0));
                        
            newBoundaries(1) = xiBoundariesAndTPs(1);
            totalBoundaryPoints = 1;
            for iInteriorPoint = 2:(length(xiBoundariesAndTPs)-1)
                % decision tree to see if we keep this interior point
                if thesign(iInteriorPoint-1) > 0
                    % positive to left...
                    if iInteriorPoint == length(xiBoundariesAndTPs)-1 &&  LxiRegion(iInteriorPoint) > L_osc/2
                        %...and bottom boundary, with non-penetrating region
                        % so make the positive (oscillatory) region bigger
                        newsigns(totalBoundaryPoints) = 1;
                        totalBoundaryPoints = totalBoundaryPoints + 1;
                        newBoundaries(totalBoundaryPoints) = xiBoundariesAndTPs(iInteriorPoint) - L_osc/2;     
                    elseif iInteriorPoint < length(xiBoundariesAndTPs)-1 &&  LxiRegion(iInteriorPoint) > L_osc
                        %...and interior boundary, with non-penetrating region
                        newsigns(totalBoundaryPoints) = 1;
                        totalBoundaryPoints = totalBoundaryPoints + 1;
                        newBoundaries(totalBoundaryPoints) = xiBoundariesAndTPs(iInteriorPoint) - L_osc/2;
                    end
                elseif thesign(iInteriorPoint-1) < 0
                    % negative to the left...
                    if iInteriorPoint == 2 && LxiRegion(iInteriorPoint-1) > L_osc/2
                        %...and top boundary, with non-penetrating region
                        newsigns(totalBoundaryPoints) = -1;
                        totalBoundaryPoints = totalBoundaryPoints + 1;
                        newBoundaries(totalBoundaryPoints) = xiBoundariesAndTPs(iInteriorPoint) + L_osc/2;
                    elseif iInteriorPoint == 2 && LxiRegion(iInteriorPoint-1) > L_osc
                        %...and interior boundary, with non-penetrating region
                        newsigns(totalBoundaryPoints) = -1;
                        totalBoundaryPoints = totalBoundaryPoints + 1;
                        newBoundaries(totalBoundaryPoints) = xiBoundariesAndTPs(iInteriorPoint) + L_osc/2;
                    end
                end
            end
            if totalBoundaryPoints>1 % the last region will always have the opposite sign as the previous
                newsigns(totalBoundaryPoints) = -1*newsigns(totalBoundaryPoints-1);
            else
                newsigns(totalBoundaryPoints) = 1;
            end
            totalBoundaryPoints = totalBoundaryPoints + 1;
            newBoundaries(totalBoundaryPoints) = xiBoundariesAndTPs(end);

            
            newsigns = newsigns(1:(totalBoundaryPoints-1));
            newBoundaries = newBoundaries(1:totalBoundaryPoints);
        end
        
        function [newBoundaries, newSigns] = RemoveRegionAtIndex(self, oldBoundaries, oldSigns, index)
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
        
        function nEVPPoints = DistributePointsInRegions(self, nTotalPoints, boundaries, thesign)
            L = abs(diff(boundaries));
            totalEquations = length(boundaries)-1;
            
            % This algorithm distributes the points/polynomials amongst the
            % different regions/equations
            if totalEquations > 1
                ratio = 1/6;
                
                % Never let a region get more points than it would have,
                % and cap that ratio.
                ratio = min( sum(L(thesign<0))/sum(L), ratio );
                
                nNegativePoints = floor(ratio*nTotalPoints);
                nPositivePoints = nTotalPoints - nNegativePoints;
                nEVPPoints = zeros(totalEquations,1);
                
                indices = thesign<0;
                if ~isempty(indices)
                    relativeSize = sqrt( L(indices)./min(L(indices)));
                    nBase = nNegativePoints/sum( relativeSize );
                    nEVPPoints(indices) = floor( relativeSize * nBase );
                    nEVPPoints(indices(end)) = nEVPPoints(indices(end)) + nNegativePoints-sum(nEVPPoints(indices));
                end
                
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
        
        function nEVPPoints = DistributePointsInRegionsWithMinimum(self, nTotalPoints, boundaries, thesign)
            L = abs(diff(boundaries));
            totalEquations = length(boundaries)-1;
            minPoints = 6;
            
            % This algorithm distributes the points/polynomials amongst the
            % different regions/equations
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
        
        function SetupCoupledEquationsAtBoundaries( self, boundaries, nEVPPoints )
            self.xiBoundaries = boundaries;
            self.zBoundaries = interp1(self.xi_zLobatto,self.zLobatto,self.xiBoundaries,'spline');
            self.Lxi = abs(diff(boundaries));
            self.nEquations = length(self.xiBoundaries)-1;
            
            % A boundary point is repeated at the start of each EVP
            boundaryIndicesStart = cumsum( [1; nEVPPoints(1:end-1)] );
            boundaryIndicesEnd = boundaryIndicesStart + nEVPPoints-1;
            
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
            self.xiLobatto = zeros(nGridPoints,1);
            self.Int_xCheb = zeros(self.nEVP,1);
            self.T_xLobatto = zeros(self.nEVP,self.nEVP);
            self.Tx_xLobatto = zeros(self.nEVP,self.nEVP);
            self.Txx_xLobatto = zeros(self.nEVP,self.nEVP);
            for i=1:self.nEquations
                n = length(self.eqIndices{i});
                m = length(self.polyIndices{i});
                xLobatto = (self.Lxi(i)/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + min(self.xiBoundaries(i+1),self.xiBoundaries(i));
                self.xiLobatto(self.xiIndices{i}) = xLobatto;
                
                [T,Tx,Txx] = InternalModesSpectral.ChebyshevPolynomialsOnGrid( xLobatto, m );
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
                    
                    self.T_xCheb_zOut_fromIndices{index} = self.polyIndices{i}(1:length(self.xiIndices{i}));
                    self.T_xCheb_zOut_toIndices{index} = toIndices;
                    self.T_xCheb_zOut_Transforms{index} = InternalModesSpectral.ChebyshevTransformForGrid(self.xiLobatto(self.xiIndices{i}), self.xiOut(toIndices));
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
