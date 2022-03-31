classdef InternalModesFiniteDifference < InternalModesBase
    % This class uses finite differencing of arbitrary order to compute the
    % internal wave modes. See InternalModesBase for basic usage
    % information.
    %
    %   The class takes the name/value pair 'orderOfAccuracy' (default
    %   value of 4) in order to set the order of accuracy of the finite
    %   differencing matrix. The matrix is constructed using the weights
    %   algorithm described by Bengt Fornberg in 'Calculation of weight in
    %   finite difference formulas', SIAM review, 1998.
    %
    %   Setting the orderOfAccuracy does tended to improve the quality of
    %   the solution, but does tend to have strange effects when the order
    %   gets high relative to the number of grid points.
    %
    %   If you request a different output grid than input grid, the
    %   solutions are mapped to the output grid with spline interpolation.
    %
    %   See also INTERNALMODES, INTERNALMODESSPECTRAL,
    %   INTERNALMODESDENSITYSPECTRAL, INTERNALMODESWKBSPECTRAL, and
    %   INTERNALMODESBASE.
    %
    %
    %   Jeffrey J. Early
    %   jeffrey@jeffreyearly.com
    %
    %   March 14th, 2017        Version 1.0
    
    properties (Access = public)
        rho  % Density on the z grid.
        N2 % Buoyancy frequency on the z grid, $N^2 = -\frac{g}{\rho(0)} \frac{\partial \rho}{\partial z}$.
        orderOfAccuracy = 4 % Order of accuracy of the finite difference matrices.
    end
    
    properties (Dependent)
        rho_z % First derivative of density on the z grid.
        rho_zz % Second derivative of density on the z grid.
    end
    
    properties (Access = public)
        n                   % length of z_diff
        z_diff              % the z-grid used for differentiation
        rho_z_diff          % rho on the z_diff grid
        N2_z_diff           % N2 on the z_diff grid
        Diff1               % 1st derivative matrix, w/ 1st derivative boundaries
        Diff2               % 2nd derivative matrix, w/ BCs set by upperBoundary property
        
        T_out               % *function* handle that transforms from z_diff functions to z_out functions
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = InternalModesFiniteDifference(rho, z_in, z_out, latitude, varargin)
            % Initialize with either a grid or analytical profile.
            self@InternalModesBase(rho,z_in,z_out,latitude,varargin{:});
            
            self.n = length(self.z_diff);
            self.Diff1 = InternalModesFiniteDifference.FiniteDifferenceMatrix(1, self.z_diff, 1, 1, self.orderOfAccuracy);
            self.Diff2 = InternalModesFiniteDifference.FiniteDifferenceMatrix(2, self.z_diff, 2, 2, self.orderOfAccuracy);
            self.N2_z_diff = -(self.g/self.rho0) * self.Diff1 * self.rho_z_diff;
            
            self.InitializeOutputTransformation(z_out);
            self.rho = self.T_out(self.rho_z_diff);
            self.N2 = self.T_out(self.N2_z_diff);
            
            if isempty(self.nModes) || self.nModes < 1
                self.nModes = self.n;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computation of the modes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [F,G,h,omega,varargout] = ModesAtWavenumber(self, k, varargin )
            % Return the normal modes and eigenvalue at a given wavenumber.
            
            self.gridFrequency = 0;
            
            % The eigenvalue equation is,
            % G_{zz} - K^2 G = \frac{f_0^2 -N^2}{gh_j}G
            % A = \left( \partial_{zz} - K^2*I \right)
            % B = \frac{f_0^2 - N^2}{g}
            A = self.Diff2 - k*k*eye(self.n);
            B = diag(self.f0*self.f0 - self.N2_z_diff)/self.g;
            
            [A,B] = self.ApplyBoundaryConditions(A,B);
            
            h_func = @(lambda) 1.0 ./ lambda;
            varargout = cell(size(varargin));
            [F,G,h,varargout{:}] = ModesFromGEP(self,A,B,h_func, varargin{:});
            omega = self.omegaFromK(h,k);
        end
        
        function [F,G,h,k,varargout] = ModesAtFrequency(self, omega, varargin )
            % Return the normal modes and eigenvalue at a given frequency.
            
            self.gridFrequency = omega;
            
            A = self.Diff2;
            B = -diag(self.N2_z_diff - omega*omega)/self.g;
            
            [A,B] = self.ApplyBoundaryConditions(A,B);
                        
            h_func = @(lambda) 1.0 ./ lambda;
            varargout = cell(size(varargin));
            [F,G,h,varargout{:}] = ModesFromGEP(self,A,B,h_func, varargin{:});
            k = self.kFromOmega(h,omega);
        end
        
        function [A,B] = ApplyBoundaryConditions(self,A,B)
            iSurface = size(A,1);
            iBottom = 1;
            
            switch self.lowerBoundary
                case LowerBoundary.freeSlip % G = 0
                    A(iBottom,:) = 0;
                    A(iBottom,iBottom) = 1;
                    B(iBottom,:) = 0;
                case LowerBoundary.noSlip % G_z = 0
                    D = weights(self.z_diff(iBottom),self.z_diff,1);
                    A(iBottom,:) = D(2,:);
                    B(iBottom,:) = 0;
                case LowerBoundary.none
                otherwise
                    error('Unknown boundary condition');
            end
            
            % G=0 or N^2 G_s = \frac{1}{h_j} G at the surface, depending on the BC
            switch self.upperBoundary
                case UpperBoundary.freeSurface
                    % G_z = \frac{1}{h_j} G at the surface
                    range = (iSurface-(self.orderOfAccuracy+1-1)):iSurface;
                    D = InternalModesFiniteDifference.weights( self.z_diff(iSurface), self.z_diff(range), 1 );
                    A(iSurface,:) = 0;
                    A(iSurface,range) = D(2,:);
                    B(iSurface,:) = 0;
                    B(iSurface,iSurface) = 1;
                case UpperBoundary.rigidLid
                    A(iSurface,:) = 0;
                    A(iSurface,iSurface) = 1;
                    B(iSurface,:) = 0;
                case UpperBoundary.none
                otherwise
                    error('Unknown boundary condition');
            end
        end
        
        function psi = SurfaceModesAtWavenumber(self, k)
            psi = self.BoundaryModesAtWavenumber(k,01);
        end
        
        function psi = BottomModesAtWavenumber(self, k)
            psi = self.BoundaryModesAtWavenumber(k,0);
        end
        
        function psi = BoundaryModesAtWavenumber(self, k, isSurface)
            sizeK = size(k);
            if length(sizeK) == 2 && sizeK(2) == 1
                sizeK(2) = [];
            end
            
            % f'' finite diff matrix with f' at the boundaries
            diff2 = InternalModesFiniteDifference.FiniteDifferenceMatrix(2, self.z_diff, 1, 1, self.orderOfAccuracy);
            N2z_z_diff = -(self.g/self.rho0) * diff2 * self.rho_z_diff;
            A = self.N2_z_diff .* diff2 - N2z_z_diff .* self.Diff1;
            B = - (1/(self.f0*self.f0))* (self.N2_z_diff.*self.N2_z_diff) .* eye(self.n);
            
            b = zeros(self.n,1);
            if isSurface == 1
                b(end) = 1;
            else
                b(1) = 1;
            end

            psi = zeros(length(k),self.n);
            for ii = 1:length(k)
                M = A + k(ii)*k(ii)*B;
                M(1,:) = self.f0*diff2(1,:);
                M(end,:) = self.f0*diff2(end,:);
                psi(ii,:) = M\b;
            end
            
            sizeK(end+1) = self.n;
            psi = reshape(psi,sizeK);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computed (dependent) properties
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function value = get.rho_z(self)
            value = self.Diff1 * self.rho_z_diff;
        end
        
        function value = get.rho_zz(self)
            diff2 = InternalModesFiniteDifference.FiniteDifferenceMatrix(2, self.z_diff, 2, 2, self.orderOfAccuracy);
            value = diff2 * self.rho_z_diff;
        end
    end
    
    methods (Access = protected)     
        
        function self = InitializeWithBSpline(self, rho)
           error('Not yet implemented') 
        end
        
        function self = InitializeWithGrid(self, rho, z_in)
            % Used internally by subclasses to intialize with a density function.
            %
            % Superclass calls this method upon initialization when it
            % determines that the input is given in gridded form. The goal
            % is to initialize z_diff and rho_z_diff.
            self.z_diff = z_in;
            self.rho_z_diff = rho;
        end

        function self = InitializeWithFunction(self, rho, z_min, z_max)
            % Used internally by subclasses to intialize with a density grid.
            %
            % The superclass calls this method upon initialization when it
            % determines that the input is given in functional form. The
            % goal is to initialize z_diff and rho_z_diff.
            if length(self.z) < 5
                error('You need more than 5 point output points for finite differencing to work');
            end
            
            if (min(self.z) == z_min && max(self.z) == z_max)
                self.z_diff = self.z;
                self.rho_z_diff = rho(self.z_diff);
            else
                error('Other cases not yet implemented');
                % Eventually we may want to use stretched coordinates as a
                % default
            end
        end
        
        function self = InitializeWithN2Function(self, N2, zMin, zMax)
            fprintf('Initialization from N2 has not yet been unit tested...or implemented for that matter.');
            % Note that there will be a grid mismatch here---so we need to
            % do something clever...
%             self.z_diff = z_in;
%             self.rho_z_diff = -(self.rho0/self.g)*N2;
        end
    end
    
    methods (Access = public)   
        function self = InitializeOutputTransformation(self, z_out)
            % After the input variables have been initialized, this is used to
            % initialize the output transformation, T_out(f).            
            if isequal(self.z_diff,z_out)
                self.T_out = @(f_in) real(f_in);
            else % want to interpolate onto the output grid
                self.T_out = @(f_in) interp1(self.z_diff,real(f_in),z_out);
            end
        end
        

        function [F,G,h,varargout] = ModesFromGEP(self,A,B,h_func, varargin)
            % Take matrices A and B from the generalized eigenvalue problem
            % (GEP) and returns F,G,h. The h_func parameter is a function that
            % returns the eigendepth, h, given eigenvalue lambda from the GEP.
            [V,D] = eig( A, B );
            
            [h, permutation] = sort(real(h_func(diag(D))),'descend');
            G = V(:,permutation);
            
            F = zeros(self.n,self.n);
            for j=1:self.n
                F(:,j) = h(j) * self.Diff1 * G(:,j);
            end
            
            if isempty(varargin)
                [F_norm,G_norm] = self.NormalizeModes(F,G,self.N2_z_diff, self.z_diff);
            else
                varargout = cell(size(varargin));
                [F_norm,G_norm,varargout{:}] = self.NormalizeModes(F,G,self.N2_z_diff, self.z_diff, varargin{:});
            end
            
            F = zeros(length(self.z),self.nModes);
            G = zeros(length(self.z),self.nModes);
            for iMode=1:self.nModes
                F(:,iMode) = self.T_out(F_norm(:,iMode));
                G(:,iMode) = self.T_out(G_norm(:,iMode));
            end
            h = reshape(h(1:self.nModes),1,[]);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Generical function to normalize
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [F,G,varargout] = NormalizeModes(self,F,G,N2,z,varargin)
            % This method normalizes the modes F,G using trapezoidal
            % integration on the given z grid. At the moment, this is only
            % used by the finite differencing algorithm, as the spectral
            % methods can use a superior (more accurate) technique of
            % directly integrating the polynomials.
            if z(2)-z(1) > 0
                direction = 'last';
            else
                direction = 'first';
            end
            
            varargout = cell(size(varargin));
            for iArg=1:length(varargin)
                varargout{iArg} = zeros(1,length(G(1,:)));
            end
            
            [maxIndexZ] = find(N2-self.gridFrequency*self.gridFrequency>0,1,direction);  
            for j=1:length(G(1,:))
                switch self.normalization
                    case Normalization.uMax
                        A = max( abs(F(:,j)) );
                        G(:,j) = G(:,j) / A;
                        F(:,j) = F(:,j) / A;
                    case Normalization.wMax
                        A = max( abs(G(:,j)) );
                        G(:,j) = G(:,j) / A;
                        F(:,j) = F(:,j) / A;
                    case Normalization.kConstant
                        if z(2)-z(1) > 0
                            G20 = G(end,j)^2;
                        else
                            G20 = G(1,j)^2;
                        end
                        A = abs(G20 + trapz( z, (1/self.g) * (N2 - self.f0*self.f0) .* G(:,j) .^ 2));
                        G(:,j) = G(:,j) / sqrt(A);
                        F(:,j) = F(:,j) / sqrt(A);
                    case Normalization.omegaConstant
                        A = abs(trapz( z, (1/abs(z(end)-z(1))) .* F(:,j) .^ 2));
                        G(:,j) = G(:,j) / sqrt(A);
                        F(:,j) = F(:,j) / sqrt(A);
                end
                
                if F(maxIndexZ,j)< 0
                    F(:,j) = -F(:,j);
                    G(:,j) = -G(:,j);
                end
                
                for iArg=1:length(varargin)
                    if ( strcmp(varargin{iArg}, 'F2') )
                        varargout{iArg}(j) = abs(trapz( z, F(:,j) .^ 2));
                    elseif ( strcmp(varargin{iArg}, 'G2') )
                        varargout{iArg}(j) = abs(trapz(z, G(:,j).^2));
                    elseif ( strcmp(varargin{iArg}, 'N2G2') )
                        varargout{iArg}(j) = abs(trapz(z, N2.* (G(:,j).^2)));
                    elseif  ( strcmp(varargin{iArg}, 'uMax') )
                        B = max( abs(F(:,j)) );
                        varargout{iArg}(j) = abs(1/B);
                    elseif  ( strcmp(varargin{iArg}, 'wMax') )
                        B = max( abs(G(:,j)) );
                        varargout{iArg}(j) = abs(1/B);
                    elseif ( strcmp(varargin{iArg}, 'kConstant') )
                        if z(2)-z(1) > 0
                            G20 = G(end,j)^2;
                        else
                            G20 = G(1,j)^2;
                        end
                        B = abs(G20 + trapz( z, (1/self.g) * (N2 - self.f0*self.f0) .* G(:,j) .^ 2));
                        varargout{iArg}(j) = sqrt(abs(1/B));
                    elseif ( strcmp(varargin{iArg}, 'omegaConstant') )
                        B = abs(trapz( z, (1/abs(z(end)-z(1))) .* F(:,j) .^ 2));
                        varargout{iArg}(j) = sqrt(abs(1/B));
                    else
                        error('Invalid option. You may request F2, G2, N2G2');
                    end
                end
            end
        end
    end
    
    methods (Static)
        function D = FiniteDifferenceMatrix(numDerivs, x, leftBCDerivs, rightBCDerivs, orderOfAccuracy)
            % Creates a finite difference matrix of aribtrary accuracy, on an arbitrary
            % grid. Left and right boundary conditions are specified as their order of
            % derivative.
            %
            % numDerivs ? the number of derivatives
            % x ? the grid
            % leftBCDerivs ? derivatives for the left boundary condition.
            % rightBCDerivs ? derivatives for the right boundary condition.
            % orderOfAccuracy ? minimum order of accuracy required
            %
            % Jeffrey J. Early, 2015
            
            n = length(x);
            D = zeros(n,n);
            
            % left boundary condition
            range = 1:(orderOfAccuracy+leftBCDerivs); % not +1 because we're computing inclusive
            c = InternalModesFiniteDifference.weights( x(1), x(range), leftBCDerivs );
            D(1,range) = c(leftBCDerivs+1,:);
            
            % central derivatives, including possible weird end points
            centralBandwidth = ceil(numDerivs/2)+ceil(orderOfAccuracy/2)-1;
            for i=2:(n-1)
                rangeLength = 2*centralBandwidth; % not +1 because we're computing inclusive
                startIndex = max(i-centralBandwidth, 1);
                endIndex = startIndex+rangeLength;
                if (endIndex > n)
                    endIndex = n;
                    startIndex = endIndex-rangeLength;
                end
                range = startIndex:endIndex;
                c = InternalModesFiniteDifference.weights( x(i), x(range), numDerivs );
                D(i,range) = c(numDerivs+1,:);
            end
            
            % right boundary condition
            range = (n-(orderOfAccuracy+rightBCDerivs-1)):n; % not +1 because we're computing inclusive
            c = InternalModesFiniteDifference.weights( x(n), x(range), rightBCDerivs );
            D(n,range) = c(rightBCDerivs+1,:);
        end
                
        function c = weights(z,x,m)
            % Calculates FD weights. The parameters are:
            %  z   location where approximations are to be accurate,
            %  x   vector with x-coordinates for grid points,
            %  m   highest derivative that we want to find weights for
            %  c   array size m+1,lentgh(x) containing (as output) in
            %      successive rows the weights for derivatives 0,1,...,m.
            %
            % Taken from Bengt Fornberg
            %
            n=length(x); c=zeros(m+1,n); c1=1; c4=x(1)-z; c(1,1)=1;
            for i=2:n
                mn=min(i,m+1); c2=1; c5=c4; c4=x(i)-z;
                for j=1:i-1
                    c3=x(i)-x(j);  c2=c2*c3;
                    if j==i-1
                        c(2:mn,i)=c1*((1:mn-1)'.*c(1:mn-1,i-1)-c5*c(2:mn,i-1))/c2;
                        c(1,i)=-c1*c5*c(1,i-1)/c2;
                    end
                    c(2:mn,j)=(c4*c(2:mn,j)-(1:mn-1)'.*c(1:mn-1,j))/c3;
                    c(1,j)=c4*c(1,j)/c3;
                end
                c1=c2;
            end
            
        end
    end
    
end

