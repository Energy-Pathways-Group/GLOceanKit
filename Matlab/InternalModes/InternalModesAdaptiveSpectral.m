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
            self.N2_zLobatto = InternalModesSpectral.ifct(self.N2_zCheb);
            N_zLobatto = sqrt(self.N2_zLobatto);
            self.N_zCheb = InternalModesSpectral.fct(N_zLobatto);
            self.xi_zLobatto = cumtrapz(self.zLobatto,N_zLobatto);
            
            % Now we need z on the \xi grid
            self.z_xiLobatto = interp1(self.xi_zLobatto, self.zLobatto, self.xiLobatto, 'spline');
            
            % and z_out on the \xi grid
            self.xiOut = interp1(self.zLobatto, self.xi_zLobatto, self.z, 'spline');
        end
        
        function [zBoundary, thesign] = FindTurningPointBoundariesAtFrequency(self, omega)
            % This function returns not just the turning points, but also
            % the top and bottom boundary locations in z.
            N2Omega2_zLobatto = self.N2_zLobatto - omega*omega;
            a = N2Omega2_zLobatto; a(a>=0) = 1; a(a<0) = 0;
            turningIndices = find(diff(a)~=0);
            nTP = length(turningIndices);
            zTP = zeros(nTP,1);
            for i=1:nTP
                fun = @(z) interp1(self.zLobatto,N2Omega2_zLobatto,z,'spline');
                zTP(i) = fzero(fun,self.zLobatto(turningIndices(i)));
            end
            zBoundary = [self.zLobatto(1); zTP; self.zLobatto(end)];
            
            % what's the sign on the EVP in these regions? In this case,
            % positive indicates oscillatory, negative exponential decay
            midZ = zBoundary(1:end-1) + diff(zBoundary)/2;
            thesign = sign( interp1(self.zLobatto,N2Omega2_zLobatto,midZ,'spline') );
            % remove any boundaries of zero length
            for index = reshape(find( thesign == 0),1,[])
                thesign(index) = [];
                zBoundary(index) = [];
            end
        end
        
        function CreateGridForFrequency(self,omega)
            if self.gridFrequency == omega
                return
            end
            
            self.gridFrequency = omega;
                        
            [zBoundariesAndTPs, thesign] = self.FindTurningPointBoundariesAtFrequency(self.gridFrequency);
            xiBoundariesAndTPs = interp1(self.zLobatto,self.xi_zLobatto,zBoundariesAndTPs,'spline');
            LxiRegion = abs(diff(xiBoundariesAndTPs));
            
            % The equation boundaries are different from the turning point
            % boundaries. We need to extend the oscillatories regions into
            % the exponential decay regions.
            L_osc = 3*sum(LxiRegion(thesign>0));
            % all else being equal, increasing this length scale decreases
            % the quality of the lowest modes, decreases the highest.
                        
            self.xiBoundaries(1) = xiBoundariesAndTPs(1);
            newsigns = [];
            for iInteriorPoint = 2:(length(xiBoundariesAndTPs)-1)
                % decision tree to see if we keep this interior point
                if thesign(iInteriorPoint-1) > 0
                    % positive to left...
                    if iInteriorPoint == length(xiBoundariesAndTPs)-1 &&  LxiRegion(iInteriorPoint) > L_osc/2
                        %...and bottom boundary, with non-penetrating region
                        self.xiBoundaries(end+1) = xiBoundariesAndTPs(iInteriorPoint) - L_osc/2;
                        newsigns(end+1) = 1;
                    elseif iInteriorPoint < length(xiBoundariesAndTPs)-1 &&  LxiRegion(iInteriorPoint) > L_osc
                        %...and interior boundary, with non-penetrating region
                        self.xiBoundaries(end+1) = xiBoundariesAndTPs(iInteriorPoint) - L_osc/2;
                        newsigns(end+1) = 1;
                    end
                elseif thesign(iInteriorPoint-1) < 0
                    % negative to the left...
                    if iInteriorPoint == 2 && LxiRegion(iInteriorPoint-1) > L_osc/2
                        %...and top boundary, with non-penetrating region
                        self.xiBoundaries(end+1) = xiBoundariesAndTPs(iInteriorPoint) + L_osc/2;
                        newsigns(end+1) = -1;
                    elseif iInteriorPoint == 2 && LxiRegion(iInteriorPoint-1) > L_osc
                        %...and interior boundary, with non-penetrating region
                        self.xiBoundaries(end+1) = xiBoundariesAndTPs(iInteriorPoint) + L_osc/2;
                        newsigns(end+1) = -1;
                    end
                end
            end
            
            % Let's hack in an extra point
%             self.xiBoundaries(end+1) = min(xiBoundariesAndTPs) + (max(xiBoundariesAndTPs)-min(xiBoundariesAndTPs))/3;
            
            self.xiBoundaries(end+1) = xiBoundariesAndTPs(end);
            self.xiBoundaries = reshape(self.xiBoundaries,[],1);
            if length(newsigns)>0
                newsigns(end+1) = -1*newsigns(end);
            else
                newsigns(end+1) = 1;
            end
            
            self.Lxi = abs(diff(self.xiBoundaries));
            self.zBoundaries = interp1(self.xi_zLobatto,self.zLobatto,self.xiBoundaries,'spline');
            
            % We will be coupling nTP+1 EVPs together. We need to
            % distribute the user requested points to each of these EVPs.
            self.nEquations = length(self.xiBoundaries)-1;
            if 1 == 0
                % For this first draft, we simply evenly distribute the points.    
                nPoints = floor(self.nEVP/self.nEquations);
                nEVPPoints = nPoints*ones(self.nEquations,1);
                nEVPPoints(end) = nEVPPoints(end) + self.nEVP - nPoints*self.nEquations; % add any extra points to the end
            elseif 1 == 1
                if self.nEquations > 1
                    ratio = 0.5;
                    nNegativePoints = floor(ratio*self.nEVP);
                    nPositivePoints = self.nEVP - nNegativePoints;
                    nEVPPoints = zeros(self.nEquations,1);
                                        
                    indices = newsigns<0;
                    if ~isempty(indices)
                        relativeSize = sqrt( self.Lxi(indices)./min(self.Lxi(indices)));
                        nBase = nNegativePoints/sum( relativeSize );
                        nEVPPoints(indices) = floor( relativeSize * nBase );
                        nEVPPoints(indices(end)) = nEVPPoints(indices(end)) + nNegativePoints-sum(nEVPPoints(indices));
                    end
                    
                    indices = newsigns>0;
                    if ~isempty(indices)
                        relativeSize = sqrt( self.Lxi(indices)./min(self.Lxi(indices)));
                        nBase = nPositivePoints/sum( relativeSize );
                        nEVPPoints(indices) = floor( relativeSize * nBase );
                        nEVPPoints(indices(end)) = nEVPPoints(indices(end)) + nPositivePoints-sum(nEVPPoints(indices));
                    end
                    
                else
                    nEVPPoints = self.nEVP;
                end
                    
            elseif 1 == 1
%                 nEVPPoints = zeros(self.nEquations,1);
                minN = 10;
                
                if sum(newsigns<0) > 0
                    minLxi = min(self.Lxi(newsigns<0));
                else
                    minLxi = min(self.Lxi);
                end
                nEVPPoints = floor(minN*sqrt(self.Lxi/minLxi));
                
%                 % force the smallest non-oscillatory section to have 12 points, and
%                 % everthing else proportionally more
%                 indices = newsigns<0;
%                 nEVPPoints(indices) = floor(minN*sqrt(self.Lxi(indices)/min(self.Lxi(indices))));
%                 
%                 % and then force the smallest oscillatory section to have
%                 % 12 points as well
%                 indices = newsigns>0;
%                 nEVPPoints(indices) = floor(minN*sqrt(self.Lxi(indices)/min(self.Lxi(indices))));
                
                % Now distribute the remaining points with some ratio of
                % non-oscillatory to (oscillatory+non-oscillatory).
                ratio = 0.1;
                indices = newsigns<0;
                remainingPoints = self.nEVP-sum(nEVPPoints);
                nNegativePoints = floor( ratio*remainingPoints);                
                
                Ltotal = sum(self.Lxi(indices));
                nEVPPoints(indices) = nEVPPoints(indices) + floor(nNegativePoints*self.Lxi(indices)./Ltotal);
                
                % give any extra points to the oscillatory section
                nPositivePoints = self.nEVP - sum(nEVPPoints);
                indices = newsigns>0;
                Ltotal = sum(self.Lxi(indices));
                lastOscIndex = find(newsigns>0,1,'last');
                nEVPPoints(indices) = nEVPPoints(indices) + floor(nPositivePoints*self.Lxi(indices)./Ltotal);
                
                % give any extra points to the last oscillatory section
                nEVPPoints(lastOscIndex) = nEVPPoints(lastOscIndex) + (self.nEVP - sum(nEVPPoints));
                
                if sum(nEVPPoints) > self.nEVP
                    error('removing turning points, you did not provide enough grid points')
                end
            elseif 1 == 2
                % Here we try giving oscillatory regions 2/3s of the points
                % 2017/09/15 ? this appears worse than the above, simpler
                % method!
                nNegativePoints = floor( (1/7)*self.nEVP );
                nPositivePoints = self.nEVP - nNegativePoints;
                nEVPPoints = zeros(self.nEquations,1);
                
                indices = newsigns<0;
                Ltotal = sum(self.Lxi(indices));
                nEVPPoints(indices) = floor(nNegativePoints*self.Lxi(indices)./Ltotal);
                
                % give any extra points to the oscillatory section
                nPositivePoints = nPositivePoints + (nNegativePoints - sum(nEVPPoints));
                
                indices = newsigns>0;
                Ltotal = sum(self.Lxi(indices));
                nEVPPoints(indices) = floor(nPositivePoints*self.Lxi(indices)./Ltotal);
                
                % give any extra points to the last oscillatory section
                nEVPPoints(indices(end)) = nEVPPoints(indices(end)) + (self.nEVP - sum(nEVPPoints));
            elseif 1 == 3
                nEVPPoints = floor(self.nEVP*self.Lxi/sum(self.Lxi));
                nEVPPoints(1) = nEVPPoints(1) + (self.nEVP - sum(nEVPPoints));
            elseif 1 == 1
                % The accuracy of the lowest mode is controlled by the
                % minimum grid points in the non-oscillatory region.
                % With exponential stratification test at 64 and 128 grid
                % points, we don't drop below 1e-2 error until we use 16
                % points in each section, and below 1e-3 at 24 grid points.
                %
                % The double exponential does great with 12/40/12 (2e-5)
                % and even better with 16/32/16, although at that
                % distribution the high modes look basically the same as
                % the standard wkb method.
                %
                % These are all with a decay scale of 3
                nEVPPoints = 16*ones(size(self.Lxi));
                remainingPoints = self.nEVP-sum(nEVPPoints);
                
                nNegativePoints = floor( (0/7)*remainingPoints);                
                indices = newsigns<0;
                Ltotal = sum(self.Lxi(indices));
                nEVPPoints(indices) = nEVPPoints(indices) + floor(nNegativePoints*self.Lxi(indices)./Ltotal);
                
                % give any extra points to the oscillatory section
                nPositivePoints = self.nEVP - sum(nEVPPoints);
                indices = newsigns>0;
                Ltotal = sum(self.Lxi(indices));
                lastOscIndex = find(newsigns>0,1,'last');
                nEVPPoints(indices) = nEVPPoints(indices) + floor(nPositivePoints*self.Lxi(indices)./Ltotal);
                
                % give any extra points to the last oscillatory section
                nEVPPoints(lastOscIndex) = nEVPPoints(lastOscIndex) + (self.nEVP - sum(nEVPPoints));
            end
            nEVPPoints
            
                        
            self.SetupCoupledEquationsAtBoundaries( self.xiBoundaries, nEVPPoints );
        end
        
        function SetupCoupledEquationsAtBoundaries( self, xiBoundaries, nEVPPoints )
            self.xiBoundaries = xiBoundaries;
            
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
