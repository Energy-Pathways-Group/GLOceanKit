classdef InternalModesSpectral < InternalModesBase
    % InternalModesSpectral uses Chebyshev polynomials on a z-grid to
    % compute the internal wave modes. See InternalModesBase for basic
    % usage information.
    %
    % This class takes the name/value pair 'nEVP' (default 513) in order to
    % set the default resolution of the polynomials used to solve the
    % eigenvalue problem (EVP), e.g.,
    %       modes = InternalModes(rho,zDomain,zOut,latitude, 'nEVP', 128);
    % Generally speaking, you will get half or more good quality modes (so
    % nEVP=513 should give you 250 or so quality modes, if the density
    % structure isn't too weird). Note that the time to solve the EVP
    % scales as the nEVP^3, so going higher resolution comes at great cost.
    % For best speed results, nEVP should be set to 2^n+1 to fully utilize
    % the FFT.
    %
    % All solutions are spectrally projected to the request ouput grid, and
    % therefore will remain high quality. The fastest output can be
    % achieved by using an Gauss-Lobatto grid spanning z.
    % 
    % The modes are normalized by integrating the Chebyshev polynomials on
    % the Lobatto grid using the exact integrals.
    %
    % A word on notation:
    %
    % Dimensions -- like z -- have grids, which may differ. zIn is the
    % dimension z, on some grid given by the user. zLobatto is the
    % diension z, on a Lobatto grid. zCheb would be used for the Chebyshev
    % polynomial representation of a function defined on z. Note that the
    % base class property 'z' is really zOut.
    %
    % Functions -- if they're analytical, they're defined on a dimension,
    % not a grid. If they're gridded, then the grid must also be specified.
    % Thus, the function rho(z) when analytically defined will be given as
    % either just rho (when obvious) or rho_z when necessary. If it's
    % gridded, then expect rho_zIn or rho_zLobatto.
    %
    % Transformations -- a transformation often a matrix (or generally just
    % a function) that transformation a function from one basis to another.
    % A common scenario here is that we need to go from zLobatto to zCheb.
    % In our cases we're usually going from one gridded dimension to
    % another. Inverse Chebyshev transformations are denote by T, or more
    % specifically T_zCheb_zLobatto.
    %
    % The underscore in the case of rho_z is certainly recognized as a
    % derivative (because of LaTex notation), so it is awkward with the
    % notion used here, where the underscore is used to separate the
    % function/transformation name from the grid. Thus, we leave rho_z and
    % rho_zz as the only two exceptions to the aforementioned notation
    % because they are the only public facing interface.
    %
    % The 'x' dimension used in the properties for this class are denoted
    % as such because the different subclasses uses these properties to
    % store values on different grids. In other words, this class uses 'x'
    % for the z dimension (depth), but the density subclass uses 'x' for
    % its density stretched coordinated.
    %
    %   See also INTERNALMODES, INTERNALMODESBASE,
    %   INTERNALMODESDENSITYSPECTRAL, INTERNALMODESWKBSPECTRAL, and
    %   INTERNALMODESFINITEDIFFERENCE.
    %
    %
    %   Jeffrey J. Early
    %   jeffrey@jeffreyearly.com
    %
    %   March 14th, 2017        Version 1.0
    
    properties (Access = public)

    end
    
    properties (Dependent)
        rho % Density on the z grid.
        rho_z % First derivative of density on the z grid.
        rho_zz % Second derivative of density on the z grid.
        N2 % Buoyancy frequency on the z grid, $N^2 = -\frac{g}{\rho(0)} \frac{\partial \rho}{\partial z}$.
    end
    

    properties %(Access = private)
        rho_function        % function handle to return rho at z
        N2_function         % function handle to return N2 at z
                
        nEVP = 0           % number of points in the eigenvalue problem
        
        % These properties are initialized with SetupEigenvalueProblem()
        % Most subclasses *will* override these initializations. The 'x' refers to the stretched coordinate being used.
        % This class uses x=z (depth), although they may have different numbers of points.
        
        % The 'x' refers to the stretched coordinate being used.
        % Once x_function has been set, all the properties listed below are
        % automatically created.
        x_function         % function handle to return 'x' at z (i.e., the stretched coordinate function)
        xLobatto           % stretched coordinate Lobatto grid nEVP points. z for this class, density or wkb for others.
        xDomain            % limits of the stretched coordinate [xMin xMax]
        z_xLobatto         % The value of z, at the xLobatto points
        xOut               % desired locations of the output in x-coordinate (deduced from z_out)
        Diff1_xCheb        % single derivative in spectral space, *function handle*
        T_xLobatto, Tx_xLobatto, Txx_xLobatto        % Chebyshev polys (and derivs) on the zLobatto
        T_xCheb_zOut        
        Int_xCheb          % Vector that multiplies Cheb coeffs, then sum for integral
        N2_xLobatto        % N2 on the z2Lobatto grid

        % Set on initialization by the subclass, these transformations are
        % applied after solving the EVP to transform back into z-space.
        hFromLambda;
        GOutFromVCheb;
        FOutFromVCheb;
        GFromVCheb;
        FFromVCheb;
        GNorm;
        FNorm;
        GeostrophicNorm
    end
    
    properties (Dependent)
        xMin
        xMax
        Lx
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = InternalModesSpectral(options)
            arguments
                options.rho = ''
                options.N2 function_handle = @disp
                options.zIn (:,1) double = []
                options.zOut (:,1) double = []
                options.latitude (1,1) double = 33
                options.rho0 (1,1) double {mustBePositive} = 1025
                options.nModes (1,1) double = 0
                options.nEVP = 512;
                options.rotationRate (1,1) double = 7.2921e-5
                options.g (1,1) double = 9.81
            end
            self@InternalModesBase(rho=options.rho,N2=options.N2,zIn=options.zIn,zOut=options.zOut,latitude=options.latitude,rho0=options.rho0,nModes=options.nModes,rotationRate=options.rotationRate,g=options.g);
            self.nEVP = options.nEVP;
            self.SetupEigenvalueProblem();            
        end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computed (dependent) properties
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function value = get.rho(self)
            value = self.rho_function(self.z);
        end
                
        function value = get.rho_z(self)
            rho_z_function = diff(self.rho_function);
            value = rho_z_function(self.z);
        end
        
        function value = get.rho_zz(self)
            rho_zz_function = diff(diff(self.rho_function));
            value = rho_zz_function(self.z);
        end
        
        function value = get.N2(self)
            value = self.N2_function(self.z);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computation of the modes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [A,B] = EigenmatricesForWavenumber(self, k )
            % The eigenvalue equation is,
            % G_{zz} - K^2 G = \frac{f_0^2 -N^2}{gh_j}G
            % A = \frac{g}{f_0^2 -N^2} \left( \partial_{zz} - K^2*I \right)
            % B = I
            T = self.T_xLobatto;
            Tzz = self.Txx_xLobatto;
            n = self.nEVP;
            
            A = (Tzz - k*k*eye(n)*T);
            B = diag((self.f0*self.f0-self.N2_xLobatto)/self.g)*T;
            
            [A,B] = self.ApplyBoundaryConditions(A,B);
        end
        
        function [A,B] = EigenmatricesForFrequency(self, omega )
            T = self.T_xLobatto;
            Tzz = self.Txx_xLobatto;
            
            A = Tzz;
            B = diag((omega*omega-self.N2_xLobatto)/self.g)*T;
            
            [A,B] = self.ApplyBoundaryConditions(A,B);
        end

        function [A,B] = EigenmatricesForGeostrophicGModes(self, k )
            T = self.T_xLobatto;
            Tz = self.Tx_xLobatto;
            Tzz = self.Txx_xLobatto;
            n = self.nEVP;

            A = Tzz;
            B = -diag((self.N2_xLobatto)/self.g)*T;

            A(n,:) = T(n,:)-T(1,:);
            B(n,:) = 0;

            if (self.g/(self.f0*self.f0))*(k*k) < 1
                A(1,:) = Tz(1,:) + (self.g/(self.f0*self.f0))*(k*k)*T(1,:);
            else
                A(1,:) = (self.f0*self.f0)/(self.g *k*k)* Tz(1,:) + T(1,:);
            end
            B(1,:) = 0;
        end

        function [A,B] = EigenmatricesForRigidLidGModes(self, k )
            T = self.T_xLobatto;
            Tz = self.Tx_xLobatto;
            Tzz = self.Txx_xLobatto;
            n = self.nEVP;

            A = Tzz;
            B = -diag((self.N2_xLobatto)/self.g)*T;

            A(n,:) = T(n,:)-T(1,:);
            B(n,:) = 0;

            if (self.g/(self.f0*self.f0))*(k*k) < 1
                A(1,:) = Tz(1,:) + (self.g/(self.f0*self.f0))*(k*k)*T(1,:);
            else
                A(1,:) = (self.f0*self.f0)/(self.g *k*k)* Tz(1,:) + T(1,:);
            end
            B(1,:) = 0;
        end

        function [A,B] = EigenmatricesForGeostrophicFModes(self, k )
            T = self.T_xLobatto;
            Tz = self.Tx_xLobatto;
            Tzz = self.Txx_xLobatto;
            n = self.nEVP;

            % N2z = diff(self.N2_function);
            % N2z_xLobatto = N2z(self.z_xLobatto);
            % 
            % A = diag((self.N2_xLobatto))*Tzz - diag(N2z_xLobatto)*Tz;
            % B = -diag((self.N2_xLobatto).^2/self.g)*T;

            dzLnN2 = diff(log(self.N2_function));
            dzLnN2_xLobatto = dzLnN2(self.z_xLobatto);

            A = Tzz - diag(dzLnN2_xLobatto)*Tz;
            B = -diag((self.N2_xLobatto)/self.g)*T;

            A(n,:) = Tz(n,:);
            B(n,:) = 0;

            alpha = (self.g/(self.f0*self.f0))*(k*k);
            if alpha < 1
                A(1,:) = alpha*Tz(1,:);
                B(1,:) = -B(1,:);
            else
                A(1,:) = Tz(1,:);
                B(1,:) = -B(1,:)/alpha;
            end
        end

        % function [A,B] = EigenmatricesForGeostrophicRigidLidModes(self, k )
        %     T = self.T_xLobatto;
        %     Tz = self.Tx_xLobatto;
        %     Tzz = self.Txx_xLobatto;
        %     n = self.nEVP;
        % 
        %     A = Tzz;
        %     B = -diag((self.N2_xLobatto)/self.g)*T;
        % 
        %     A(n,:) = T(n,:);
        %     B(n,:) = 0;
        % 
        %     c = (self.g/(self.f0*self.f0))*(k*k);
        %     if c < 1
        %         A(1,:) = c*Tz(1,:);
        %         B(1,:) = -(Tz(1,:) + c*T(1,:));
        % 
        %         A(n,:) = c*Tz(1,:);
        %         B(n,:) = -(Tz(1,:) + c*T(n,:));
        % 
        %         % A(n,:) = T(n,:);
        %         % B(n,:) = 0;
        %     else
        %         cinv = (self.f0*self.f0)/(self.g *k*k);
        %         A(1,:) = Tz(1,:);
        %         B(1,:) = -(cinv * Tz(1,:) + T(1,:));
        % 
        %         A(n,:) = Tz(1,:);
        %         B(n,:) = -(cinv * Tz(1,:) + T(n,:));
        % 
        %         % A(n,:) = T(n,:);
        %         % B(n,:) = 0;
        %     end
        % end

        function [A,B] = EigenmatricesForGeostrophicRigidLidGModes(self, eta0, etad)
            arguments
                self InternalModesSpectral
                eta0 (1,1) double = 0
                etad (1,1) double = 0
            end
            T = self.T_xLobatto;
            Tz = self.Tx_xLobatto;
            Tzz = self.Txx_xLobatto;
            n = self.nEVP;

            A = Tzz;
            B = -(diag(self.N2_xLobatto)/self.g)*T;

            % upper-boundary
            A(1,:) = Tz(1,:);
            B(1,:) = (eta0/self.g)*self.N2_xLobatto(1)*T(1,:);

            % lower-boundary
            A(n,:) = Tz(n,:);
            B(n,:) = (etad/self.g)*self.N2_xLobatto(n)*T(n,:);

            self.FOutFromVCheb = @(G_cheb,h) h * self.T_xCheb_zOut(self.Diff1_xCheb(G_cheb));
            self.FFromVCheb = @(G_cheb,h) h * InternalModesSpectral.ifct( self.Diff1_xCheb(G_cheb) );
        end

        function [A,B] = EigenmatricesForSmithVannesteModes(self, k )
            T = self.T_xLobatto;
            Tz = self.Tx_xLobatto;
            Tzz = self.Txx_xLobatto;
            n = self.nEVP;

            A = Tzz - diag(k*k*self.N2_xLobatto/(self.f0*self.f0))*T;
            B = -(diag(self.N2_xLobatto)/self.g)*T;
            
            % upper-boundary
            A(1,:) = Tz(1,:); %+T(1,:);
            B(1,:) = self.N2_xLobatto(1)*T(1,:)/self.g;

            % lower-boundary
            A(n,:) = Tz(n,:); %self.Lz*Tz(n,:)-T(n,:);
            B(n,:) = -self.N2_xLobatto(n)*T(n,:)/self.g;

            gamma = (self.g/(self.f0*self.f0))*(k*k);
            self.FOutFromVCheb = @(G_cheb,h) (1/(gamma + 1/h)) * self.T_xCheb_zOut(self.Diff1_xCheb(G_cheb));
            self.FFromVCheb = @(G_cheb,h) (1/(gamma + 1/h)) * InternalModesSpectral.ifct( self.Diff1_xCheb(G_cheb) );

            % % upper-boundary
            % A(1,:) = self.Lz*Tz(1,:)+T(1,:);
            % B(1,:) = 0*T(1,:);
            % 
            % % lower-boundary
            % A(n,:) = T(n,:); %self.Lz*Tz(n,:)-T(n,:);
            % B(n,:) = 0*T(n,:);
        end

        function [A,B] = EigenmatricesForMDAModes(self )
            T = self.T_xLobatto;
            Tz = self.Tx_xLobatto;
            Tzz = self.Txx_xLobatto;
            n = self.nEVP;

            A = Tzz;
            B = -(diag(self.N2_xLobatto)/self.g)*T;

            % upper-boundary
            A(1,:) = Tz(1,:); %-Tz(n,:);
            B(1,:) = 0 ;%1/self.Lz; %0*T(n,:);
            self.upperBoundary = UpperBoundary.mda;

            % lower-boundary
            A(n,:) = Tz(n,:); %self.Lz*Tz(n,:)-T(n,:);
            B(n,:) = 0; %1/self.Lz; %0*T(n,:);
            self.lowerBoundary = LowerBoundary.mda;

            % A(2,:) = T(2,:);
            % B(2,:) = 0 ;
            % 
            % % lower-boundary
            % A(n-1,:) = T(n-1,:);
            % B(n-1,:) = 0;
        end
        
        function [A,B] = ApplyBoundaryConditions(self,A,B)
            T = self.T_xLobatto;
            Tz = self.Tx_xLobatto;
            n = self.nEVP;
            
            switch self.lowerBoundary
                case LowerBoundary.freeSlip
                    A(n,:) = T(n,:);
                    B(n,:) = 0;
                case LowerBoundary.noSlip
                    A(n,:) = Tz(n,:);
                    B(n,:) = 0;
%                     A(n-1,:) = T(n,:);
%                     B(n-1,:) = 0;
                case LowerBoundary.buoyancyAnomaly
                    A(n,:) = T(n,:);
                    B(n,:) = 1;
                case LowerBoundary.custom
%                     A(n,:) = Tz(n,:);% + (self.N2_xLobatto(n)/self.g).*T(n,:);
%                     B(n,:) = -T(n,:);
%                     A(n,:) = -Tz(1,:)+Tz(n,:);
                    A(n,:) = Tz(n,:);
                    B(n,:) = 0;
                case LowerBoundary.none
                otherwise
                    error('Unknown boundary condition');
            end
            
            switch self.upperBoundary
                case UpperBoundary.freeSurface
                    % G_z = \frac{1}{h_j} G at the surface
                    A(1,:) = Tz(1,:);
                    B(1,:) = T(1,:);
                case UpperBoundary.rigidLid
                    A(1,:) = T(1,:);
                    B(1,:) = 0;
                case UpperBoundary.buoyancyAnomaly
                    A(1,:) = T(1,:);
                    B(1,:) = 1;
                case UpperBoundary.geostrophicFreeSurface
                    A(1,:) = Tz(1,:) + (self.g/(self.f0*self.f0))*(options.k*options.k)*T(1,:);
                    B(1,:) = 0;
                case UpperBoundary.custom
                    A(1,:) = T(1,:)-T(n,:);
                    A(1,:) = Tz(1,:);
                    B(1,:) = 0;
                case UpperBoundary.none
                otherwise
                    error('Unknown boundary condition');
            end
        end
        
        function [F,G,h] = MDAModes(self )
            self.gridFrequency = 0;

            [A,B] = self.EigenmatricesForMDAModes();

            [F,G,h] = self.ModesFromGEP(A,B);
        end

        function [F,G,h] = GeostrophicModesAtWavenumber(self, k )
            self.gridFrequency = 0;

            [A,B] = self.EigenmatricesForGeostrophicGModes(k);

            [F,G,h] = self.ModesFromGEP(A,B);
        end

        function [F,G,h] = GeostrophicFModesAtWavenumber(self, k )
            self.gridFrequency = 0;

            [A,B] = self.EigenmatricesForGeostrophicFModes(k);

            [F,G,h] = self.ModesFromGEP(A,B);
        end

        function [F,G,h] = GeostrophicRigidLidModes(self, eta0, etad )
            arguments
                self InternalModesSpectral
                eta0 (1,1) double = 0
                etad (1,1) double = 0
            end
            self.gridFrequency = 0;

            [A,B] = self.EigenmatricesForGeostrophicRigidLidGModes(eta0, etad);
            if eta0 < 0
                negativeEigenvalues = 1;
            else
                negativeEigenvalues = 0;
            end
            [F,G,h] = self.ModesFromGEP(A,B,negativeEigenvalues=negativeEigenvalues);
        end

        function [F,G,h] = GeostrophicSmithVannesteModesAtWavenumber(self, k )
            self.gridFrequency = 0;

            [A,B] = self.EigenmatricesForSmithVannesteModes(k);

            [F,G,h] = self.ModesFromGEP(A,B);
        end

        function [F,G,h,omega,varargout] = ModesAtWavenumber(self, k, varargin )
            self.gridFrequency = 0;
            
            [A,B] = self.EigenmatricesForWavenumber(k);
            
            if isempty(varargin)
                [F,G,h] = self.ModesFromGEP(A,B);
            else
                varargout = cell(size(varargin));
                [F,G,h,varargout{:}] = self.ModesFromGEP(A,B, varargin{:});
            end

            omega = self.omegaFromK(h,k);
        end
        
        function [F,G,h,k,varargout] = ModesAtFrequency(self, omega, varargin )
            self.gridFrequency = omega;
            
            [A,B] = self.EigenmatricesForFrequency(omega);
            
            if isempty(varargin)
                [F,G,h] = self.ModesFromGEP(A,B);
            else
                varargout = cell(size(varargin));
                [F,G,h,varargout{:}] = self.ModesFromGEP(A,B, varargin{:});
            end
            
            k = self.kFromOmega(h,omega);
        end 
        
        function psi = SurfaceModesAtWavenumber(self, k) 
            psi = self.BoundaryModesAtWavenumber(k,1);
        end
        
        function psi = BottomModesAtWavenumber(self, k) 
            psi = self.BoundaryModesAtWavenumber(k,0);
        end
        
        function psi = BoundaryModesAtWavenumber(self, k, isSurface)            
            % Estimate the grid resolution necessary to resolve the
            % smallest mode.
            sizeK = size(k);
            if length(sizeK) == 2 && sizeK(2) == 1
                sizeK(2) = [];
            end
            K = k(:);
            maxK = max(K(K>0));
            maxN = sqrt(max(self.N2(1),self.N2(end)));
            minDz = self.f0/(maxK*maxN)/10;
            n = max(128,2^ceil(log2((pi/2)*sqrt(self.Lz/minDz)+1))+1);
%             fprintf('Using %d grid points for the SQG mode\n',n);
            
            % Now create this new grid
            yLobatto = (self.Lz/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + self.zMin;
            T_yCheb_zOut = InternalModesSpectral.ChebyshevTransformForGrid(yLobatto, self.z);
            
            [T,Ty,Tyy] = InternalModesSpectral.ChebyshevPolynomialsOnGrid( yLobatto, length(yLobatto) );
            N2_yLobatto = self.N2_function(yLobatto);
            N2z_function = diff(self.N2_function);
            N2_z_yLobatto = N2z_function(yLobatto);
            
            A = N2_yLobatto .* Tyy - N2_z_yLobatto.*Ty;
            B = - (1/(self.f0*self.f0))* (N2_yLobatto.*N2_yLobatto) .*T;
            
            b = zeros(size(yLobatto));
            if isSurface == 1
                b(1) = 1;
            else
                b(end) = 1;
            end
            
            psi = zeros(length(k),length(self.z));
            for ii = 1:length(k)
                M = A + k(ii)*k(ii)*B;
                M(1,:) = self.f0*Ty(1,:);
                M(end,:) = self.f0*Ty(end,:);
                
                psi_cheb = M\b;
                psi(ii,:) = T_yCheb_zOut(psi_cheb);
            end
            
            sizeK(end+1) = length(self.z);
            psi = reshape(psi,sizeK);
        end
        
        function z_g = GaussQuadraturePointsForModesAtFrequency(self,nPoints,omega)
            [A,B] = self.EigenmatricesForFrequency(omega);
            z_g = self.GaussQuadraturePointsForEigenmatrices(nPoints,A,B);
        end
        function z_g = GaussQuadraturePointsForModesAtWavenumber(self,nPoints,k)
            [A,B] = self.EigenmatricesForWavenumber(k);
            z_g = self.GaussQuadraturePointsForEigenmatrices(nPoints,A,B);
        end
        function z_g = GaussQuadraturePointsForMDAModes(self,nPoints)
            [A,B] = self.EigenmatricesForMDAModes();
            z_g = self.GaussQuadraturePointsForEigenmatrices(nPoints,A,B);
        end
        
        function z_g = GaussQuadraturePointsForEigenmatrices(self,nPoints,A,B)
            % Now we just need to find the roots of the n+1 mode.
            % For constant stratification this should give back the
            % standard Fourier modes, i.e., an evenly spaced grid.
            %
            % Note that if the boundary conditions are such that G(0)=0 and
            % G(-D)=0, then those two points do not encode any information.
            % As such, only the first (nPoints-2) modes will encode any
            % useful information. So we'd expect cond(G(:,1:(nPoints-2))))
            % to be good (low), but not the next.
            if 2*nPoints < self.nEVP
               if ( any(any(isnan(A))) || any(any(isnan(B))) )
                   error('GLOceanKit:NaNInMatrix', 'EVP setup fail. Found at least one nan in matrices A and B.\n');
               end
               [V,D] = eig( A, B );
               
               [h, permutation] = sort(real(self.hFromLambda(diag(D))),'descend');
               G_cheb=V(:,permutation);
               maxModes = ceil(find(h>0,1,'last')/2);
               
               if maxModes < (nPoints+1)
                   error('GLOceanKit:NeedMorePoints', 'Returned %d valid modes (%d quadrature points requested) using nEVPs=%d.',maxModes,nPoints,self.nEVP);
               end
  
               % Could compute the roots of the F-modes, but nah.
%                F = self.Diff1_xCheb(G_cheb(:,nPoints-1));
%                roots = InternalModesSpectral.FindRootsFromChebyshevVector(F(1:end-1), self.z_xLobatto);
%                z_g = cat(1,min(self.z_xLobatto),reshape(roots,[],1),max(self.z_xLobatto));

               % depending on the boundary conditions and particular
               % problem, the nth mode might contain (n-1), (n), or (n+1)
               % zero crossings. If nPoints are request, we want to include
               % the boundaries in that number.
               if self.upperBoundary == UpperBoundary.mda
                    rootMode = nPoints-2;
               elseif self.upperBoundary == UpperBoundary.rigidLid
                   % n-th mode has n+1 zeros (including boundaries)
                   rootMode = nPoints-1;
               elseif self.upperBoundary == UpperBoundary.freeSurface && self.lowerBoundary == LowerBoundary.noSlip
                   % n-th mode has n zeros (including zero at lower
                   % boundary, and not zero at upper)
                   rootMode = nPoints-1;
               elseif self.upperBoundary == UpperBoundary.freeSurface
                   % n-th mode has n zeros (including zero at lower
                   % boundary, and not zero at upper)
                   rootMode = nPoints;
               end
               rootsVar = InternalModesSpectral.FindRootsFromChebyshevVector(G_cheb(:,rootMode), self.xDomain);

               % First we make sure the roots are within the bounds
               rootsVar(rootsVar<self.xMin) = self.xMin;
               rootsVar(rootsVar>self.xMax) = self.xMax;
               
               % Add the boundary points---if they are redundant, they will
               % get eliminated below
               rootsVar = cat(1,self.xMin,rootsVar,self.xMax);

               % Then we eliminate any repeats (it happens)
               rootsVar = unique(rootsVar,'stable');
               
               while (length(rootsVar) > nPoints)
                   rootsVar = sort(rootsVar);
                   F = InternalModesSpectral.IntegrateChebyshevVector(G_cheb(:,rootMode));
                   value = InternalModesSpectral.ValueOfFunctionAtPointOnGrid( rootsVar, self.xDomain, F );
                   dv = diff(value);
                   [~,minIndex] = min(abs(dv));
                   rootsVar(minIndex+1) = [];
               end

               if length(rootsVar) < nPoints
                   error('GLOceanKit:NeedMorePoints', 'Returned %d unique roots (requested %d). Maybe need more EVP.', length(rootsVar),nPoints);
               end

               z_g = reshape(rootsVar,[],1);          
               z_g = InternalModesSpectral.fInverseBisection(self.x_function,z_g,min(self.zDomain),max(self.zDomain),1e-12);
            else
                error('GLOceanKit:NeedMorePoints', 'You need at least twice as many nEVP as points you request');
            end
        end
        
        function value = get.xMin(self)
            value = self.xDomain(1);
        end
        
        function value = get.xMax(self)
            value = self.xDomain(2);
        end
        
        function value = get.Lx(self)
            value = self.xMax - self.xMin;
        end
        
        function set.x_function(self,s)
            self.x_function = s;
            
            self.recomputeStretchedGrid(s);
        end
        
        
    end
    
    methods (Access = protected)
        
        function self = recomputeStretchedGrid(self,s)
            self.xDomain = [s(self.zMin) s(self.zMax)];
            self.xLobatto = ((self.xMax-self.xMin)/2)*( cos(((0:self.nEVP-1)')*pi/(self.nEVP-1)) + 1) + self.xMin;
            [self.z_xLobatto, self.xOut] = InternalModesSpectral.StretchedGridFromCoordinate( s, self.xLobatto, self.zDomain, self.z);
                        
            self.Diff1_xCheb = @(v) (2/self.Lx)*InternalModesSpectral.DifferentiateChebyshevVector( v );
            [self.T_xLobatto,self.Tx_xLobatto,self.Txx_xLobatto] = InternalModesSpectral.ChebyshevPolynomialsOnGrid( self.xLobatto, length(self.xLobatto) );
            self.T_xCheb_zOut = InternalModesSpectral.ChebyshevTransformForGrid(self.xLobatto, self.xOut);
            
            % We use that \int_{-1}^1 T_n(x) dx = \frac{(-1)^n + 1}{1-n^2}
            % for all n, except n=1, where the integral is zero.
            np = (0:(self.nEVP-1))';
            self.Int_xCheb = -(1+(-1).^np)./(np.*np-1);
            self.Int_xCheb(2) = 0;
            self.Int_xCheb = self.Lx/2*self.Int_xCheb;
            
            self.N2_xLobatto = self.N2_function(self.z_xLobatto);
      
            if self.shouldShowDiagnostics == 1
                fprintf(' The eigenvalue problem will be solved with %d points.\n', length(self.xLobatto));
            end
        end
        
        function self = SetupEigenvalueProblem(self)
            % Subclasses will override this function.
            self.x_function = @(z) z;             

            self.hFromLambda = @(lambda) 1.0 ./ lambda;
            self.GOutFromVCheb = @(G_cheb,h) self.T_xCheb_zOut(G_cheb);
            self.FOutFromVCheb = @(G_cheb,h) h * self.T_xCheb_zOut(self.Diff1_xCheb(G_cheb));
            self.GFromVCheb = @(G_cheb,h) InternalModesSpectral.ifct(G_cheb);
            self.FFromVCheb = @(G_cheb,h) h * InternalModesSpectral.ifct( self.Diff1_xCheb(G_cheb) );
            self.GNorm = @(Gj) abs(Gj(1)*Gj(1) + sum(self.Int_xCheb .*InternalModesSpectral.fct((1/self.g) * (self.N2_xLobatto - self.f0*self.f0) .* Gj .^ 2)));
            self.GeostrophicNorm = @(Gj) abs(sum(self.Int_xCheb .*InternalModesSpectral.fct((1/self.g) * self.N2_xLobatto .* Gj .^ 2)));
            self.FNorm = @(Fj) abs(sum(self.Int_xCheb .*InternalModesSpectral.fct((1/self.Lz) * Fj.^ 2)));
        end
        
        function self = InitializeWithBSpline(self, rho)
           self.validateInitialModeAndEVPSettings();
           
           if self.requiresMonotonicDensity == 1
               z = BSpline.PointsOfSupport(rho.t_knot,rho.K);
               self.rho_function = ConstrainedSpline(z,rho.ValueAtPoints(z),rho.K,rho.t_knot,NormalDistribution(1),struct('global',ShapeConstraint.monotonicDecreasing));
           else
               self.rho_function = rho;
           end
           
           
           self.N2_function = (-self.g/self.rho0)*diff(self.rho_function);
        end
        
        function self = InitializeWithGrid(self, rho, zIn)
            self.validateInitialModeAndEVPSettings();

            K = 6; 
            if self.requiresMonotonicDensity == 1
                z_knot = InterpolatingSpline.KnotPointsForPoints(zIn,K,1);
                rho_interpolant = ConstrainedSpline(zIn,rho,K,z_knot,NormalDistribution(1),struct('global',ShapeConstraint.monotonicDecreasing));
                if self.shouldShowDiagnostics == 1
                    fprintf('Creating a %d-order monotonic spline from the %d points.\n', K, length(rho));
                end
            else
                rho_interpolant = InterpolatingSpline(zIn,rho,'K',K);
                if self.shouldShowDiagnostics == 1
                    fprintf('Creating a %d-order spline from the %d points.\n', K, length(rho));
                end
            end
            
            self.rho_function = rho_interpolant;
            self.N2_function = (-self.g/self.rho0)*diff(self.rho_function);
        end
        
        function self = InitializeWithFunction(self, rho, zMin, zMax)
            self.validateInitialModeAndEVPSettings();
            
            if ~exist('chebfun','class')
               error('The package chebfun is required when initializing with a function.')
            end
            
            if isa(rho,'chebfun')
                self.rho_function = rho;
            else
                self.rho_function = chebfun(rho,[zMin zMax]);
            end
            self.N2_function = -(self.g/self.rho0)*diff(self.rho_function);
            
            if self.requiresMonotonicDensity == 1
                % Taken from inv as part of chebfun.
                f = self.rho_function;
                fp = diff(f);
                tPoints = roots(fp);
                tol = eps;
                if ( ~isempty(tPoints) )
                    endtest = zeros(length(tPoints), 1);
                    for k = 1:length(tPoints)
                        endtest(k) = min(abs(tPoints(k) - f.domain));
                    end
                    if ( any(endtest > 100*abs(feval(f, tPoints))*tol) )
                        error('Density must be monotonic to use this class. The function you provided is not monotonic.');
                    end
                end
            end
            
            if self.shouldShowDiagnostics == 1
                fprintf('Projected the function onto %d Chebyshev polynomials\n', length(self.rho_function));
            end
        end
        
        function self = InitializeWithN2Function(self, N2, zMin, zMax)
            self.validateInitialModeAndEVPSettings();
            
            if ~exist('chebfun','class')
                error('The package chebfun is required when initializing with a function.')
            end
            
            N2_func = chebfun(N2,[zMin,zMax],'splitting','on');
            rho_func = -(self.rho0/self.g)*cumsum(N2_func);
            rho_func = rho_func - rho_func(zMax) + self.rho0;
            
            self.InitializeWithFunction(rho_func,zMin,zMax);
        end
        
                   
        function self = validateInitialModeAndEVPSettings(self)
            % The user requested that the eigenvalue problem be solved on a
            % grid of particular length
            if self.nEVP > 0
                if self.nModes > self.nEVP
                    self.nEVP = self.nModes;
                end
            else
                self.nEVP = 513; % 2^n + 1 for a fast Chebyshev transform
            end            
        end
                        
        
        % Take matrices A and B from the generalized eigenvalue problem
        % (GEP) and returns F,G,h. The last seven arguments are all
        % function handles that do as they say.
        function [F,G,h,varargout] = ModesFromGEP(self,A,B,varargin,options)
            arguments
                self InternalModesSpectral
                A (:,:) double
                B (:,:) double
            end
            arguments (Repeating)
                varargin
            end
            arguments
                options.negativeEigenvalues = 0
            end
            if ( any(any(isnan(A))) || any(any(isnan(B))) )
                error('EVP setup fail. Found at least one nan in matrices A and B.\n');
            end
            [V,D] = eig( A, B );

            % The following might be better, as it captures the the
            % barotopic mode.
            % d = diag(D);
            % [d, permutation] = sort(real(d),'ascend');
            % V_cheb=V(:,permutation);
            % if options.negativeEigenvalues > 0
            %     negIndices = find(d<0,options.negativeEigenvalues,'last');
            % else
            %     negIndices = [];
            % end
            % if isempty(negIndices)
            %     minIndex = find(d>=0,1,'first');
            %     if isempty(minIndex)
            %         fprintf('No usable modes found! Try with higher resolution.\n');
            %         return;
            %     end
            % else
            %     minIndex = min(negIndices);
            % end
            % V_cheb = V_cheb(:,minIndex:end);
            % h = self.hFromLambda(d(minIndex:end));

            [h, permutation] = sort(real(self.hFromLambda(diag(D))),'descend');
            V_cheb=V(:,permutation);
            if options.negativeEigenvalues > 0
                negIndices = find(h<0,options.negativeEigenvalues,'first');
                permutation = cat(1,negIndices,setdiff((1:length(h))',negIndices));
                h = h(permutation);
                V_cheb=V_cheb(:,permutation);
            end
            
            if self.nModes == 0
                maxModes = ceil(find(h>0,1,'last')/2); % Have to do ceil, not floor, or we lose the barotropic mode.
                if maxModes == 0
                    fprintf('No usable modes found! Try with higher resolution.\n');
                    return;
                end
            else
                maxModes = self.nModes;
            end
            
%             FzOutFromGCheb = @(G_cheb,h) h * self.T_xCheb_zOut(self.Diff1_xCheb(self.Diff1_xCheb(G_cheb)));
%             Fz = zeros(length(self.z),maxModes);

            F = zeros(length(self.z),maxModes);
            G = zeros(length(self.z),maxModes);
            h = reshape(h(1:maxModes),1,[]);
            
            varargout = cell(size(varargin));
            for iArg=1:length(varargin)
                varargout{iArg} = zeros(1,maxModes);
            end
            
            % This still need to be optimized to *not* do the transforms
            % twice, when the EVP grid is the same as the output grid.
            [maxIndexZ] = find(self.N2_xLobatto-self.gridFrequency*self.gridFrequency>0,1,'first');
            if maxIndexZ > 1 % grab a point just above the turning point, which should have the right sign.
                maxIndexZ = maxIndexZ-1;
            elseif isempty(maxIndexZ)
                maxIndexZ = 1;
            end
            for j=1:maxModes
                Fj = self.FFromVCheb(V_cheb(:,j),h(j));
                Gj = self.GFromVCheb(V_cheb(:,j),h(j));
                switch self.normalization
                    case Normalization.uMax
                        A = max( abs( Fj ));
                    case Normalization.wMax
                        A = max( abs( Gj ) );
                    case Normalization.kConstant
                        A = sqrt(self.GNorm( Gj ));
                    case Normalization.omegaConstant
                        A = sqrt(self.FNorm( Fj ));
                    case Normalization.geostrophic
                        A = sqrt(self.GeostrophicNorm( Gj ));
                end
                if Fj(maxIndexZ) < 0
                    A = -A;
                end
                
                G(:,j) = self.GOutFromVCheb(V_cheb(:,j),h(j))/A;
                F(:,j) = self.FOutFromVCheb(V_cheb(:,j),h(j))/A;
%                 Fz(:,j) = FzOutFromGCheb(G_cheb(:,j),h(j))/A;
                % K-constant norm: G(0)^2 + \frac{1}{g} \int_{-D}^0 (N^2 -
                % f_0^2)
                for iArg=1:length(varargin)
                    if ( strcmp(varargin{iArg}, 'F2') )
                        varargout{iArg}(j) = self.Lz*self.FNorm( Fj/A );
                    elseif ( strcmp(varargin{iArg}, 'G2') )
                        varargout{iArg}(j) = self.Lz*self.FNorm( Gj/A );
                    elseif ( strcmp(varargin{iArg}, 'N2G2') )
                        varargout{iArg}(j) = self.g*(self.GNorm( Gj/A )-Gj(1)*Gj(1)) + self.f0*self.f0*self.Lz*self.FNorm( Gj/A ); % this is being clever, but should give \int N2*G2 dz
                    elseif  ( strcmp(varargin{iArg}, 'uMax') )
                        B = max( abs( Fj ));
                        varargout{iArg}(j) = abs(A/B);
                    elseif  ( strcmp(varargin{iArg}, 'wMax') )
                        B = max( abs( Gj ) );
                        varargout{iArg}(j) = abs(A/B);
                    elseif ( strcmp(varargin{iArg}, 'kConstant') )
                        B = sqrt(self.GNorm( Gj ));
                        varargout{iArg}(j) = abs(A/B);
                    elseif ( strcmp(varargin{iArg}, 'omegaConstant') )
                        B = sqrt(self.FNorm( Fj ));
                        varargout{iArg}(j) = abs(A/B);
                    elseif ( strcmp(varargin{iArg}, 'geostrophicNorm') )
                        B = sqrt(self.GeostrophicNorm( Gj ));
                        varargout{iArg}(j) = abs(A/B);
                    elseif ( strcmp(varargin{iArg}, 'int_N2_G_dz/g') )
                        varargout{iArg}(j) = sum(self.Int_xCheb .*InternalModesSpectral.fct((1/self.g) * self.N2_xLobatto .* (Gj/A)));
                    elseif ( strcmp(varargin{iArg}, 'int_F_dz') )
                        varargout{iArg}(j) = sum(self.Int_xCheb .*InternalModesSpectral.fct(Fj/A));
                    elseif ( strcmp(varargin{iArg}, 'int_G_dz') )
                        varargout{iArg}(j) = sum(self.Int_xCheb .*InternalModesSpectral.fct(Gj/A));
                    else
                        error('Invalid option. You may request F2, G2, N2G2');
                    end
                end
            end
            

        end
        
%         function zTPs = FindTurningPointsAtFrequency(self, omega)
%             f_cheb = self.N2_zCheb;
%             f_cheb(1) = f_cheb(1) - omega*omega;
%             zTPs= InternalModesSpectral.FindRootsFromChebyshevVector(f_cheb,self.zLobatto);
%         end
    end
    
    methods (Static)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Methods to find the turning points (used for normalization)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [zBoundariesAndTPs, thesign, boundaryIndices] = FindTurningPointBoundariesAtFrequency(N2, z, omega)
            % This function returns not just the turning points, but also
            % the top and bottom boundary locations in z. The boundary
            % indices are the index to the point just *above* the turning
            % point.
            N2Omega2 = N2 - omega*omega;
            a = N2Omega2; a(a>=0) = 1; a(a<0) = 0;
            turningIndices = find(diff(a)~=0);
            nTP = length(turningIndices);
            zTP = zeros(nTP,1);
            for i=1:nTP
                fun = @(depth) interp1(z,N2Omega2,depth,'spline');
                z0 = [z(turningIndices(i)+1) z(turningIndices(i))];
                zTP(i) = fzero(fun,z0);
            end
            boundaryIndices = [1; turningIndices; length(z)];
            zBoundariesAndTPs = [z(1); zTP; z(end)];
            
            % what's the sign on the EVP in these regions? In this case,
            % positive indicates oscillatory, negative exponential decay
            midZ = zBoundariesAndTPs(1:end-1) + diff(zBoundariesAndTPs)/2;
            thesign = sign( interp1(z,N2Omega2,midZ,'spline') );
            % remove any boundaries of zero length
            for index = reshape(find( thesign == 0),1,[])
                thesign(index) = [];
                zBoundariesAndTPs(index) = [];
                boundaryIndices(index) = [];
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Convert to a stretched grid
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [z_xLobatto, xOut] = StretchedGridFromCoordinate( x, xLobatto, zDomain, zOut)
            % x is a function handle, that maps from z.
            % xLobatto is the grid on x
            % zLobatto is the Lobatto grid on z
            % zOut is the output grid
            
            maxZ = max(zDomain);
            minZ = min(zDomain);

            z_xLobatto = InternalModesSpectral.fInverseBisection(x,xLobatto,min(zDomain),max(zDomain),1e-12);
            z_xLobatto(z_xLobatto>maxZ) = maxZ;
            z_xLobatto(z_xLobatto<minZ) = minZ;
            
            z_xLobatto(1) = maxZ;
            z_xLobatto(end) = minZ;
            
                        
            xOut = x(zOut);
            xOut(xOut>max(xLobatto)) = max(xLobatto);
            xOut(xOut<min(xLobatto)) = min(xLobatto);
        end
        
        function [flag, dTotalVariation, rho_zCheb, rho_zLobatto, rhoz_zCheb, rhoz_zLobatto] = CheckIfReasonablyMonotonic(zLobatto, rho_zCheb, rho_zLobatto, rhoz_zCheb, rhoz_zLobatto)
            % We want to know if the density function is decreasing as z
            % increases. If it's not, are the discrepencies small enough
            % that we can just force them?
            %
            % flag is 0 if the function is not reasonably monotonic, 1 if
            % it is, and 2 if we were able to to coerce it to be, without
            % too much change
            
            flag = 0;
            dTotalVariation = 0;
            if any(rhoz_zLobatto > 0)
                % record the density at the bottom, and the total variation
                % in density
                rho_top = min(rho_zLobatto(1),rho_zLobatto(end));
                rho_bottom = max(rho_zLobatto(1),rho_zLobatto(end));
                dRho = rho_bottom-rho_top;
                Lz = abs(zLobatto(1)-zLobatto(end));
                
                % Now zero out the all the regions where there are density
                % inversions, project onto cheb basis, the integrate to get
                % a new density function, but keep the density at the
                % surface the same (because we use rho0 elsewhere).
                rhoz_zLobatto(rhoz_zLobatto >= 0) = max(rhoz_zLobatto(rhoz_zLobatto<0));
                rhoz_zCheb = InternalModesSpectral.fct(rhoz_zLobatto);
                rho_zCheb = (Lz/2)*InternalModesSpectral.IntegrateChebyshevVector(rhoz_zCheb);
                rho_zCheb(end) = [];
                rho_zLobatto = InternalModesSpectral.ifct(rho_zCheb);
                delta = - min(rho_zLobatto) + rho_top;
                rho_zLobatto = rho_zLobatto + delta;
                rho_zCheb(1) = rho_zCheb(1) + delta;
                
                % Re-derive all the properties
%                 new_rhoz_zLobatto = InternalModesSpectral.ifct((2/Lz)*InternalModesSpectral.DifferentiateChebyshevVector(rho_zCheb));
                
                
                
                dRho_new = abs(rho_zLobatto(end)-rho_zLobatto(1));
                dTotalVariation = (dRho_new-dRho)/dRho;
                
                if dTotalVariation < 1e-2
                    flag = 1;
                else
                    flag = 2;
                end
            end
            
            
        end
        
        function y = fInverseBisection(f, x, yMin,yMax, tol)
            %FINVERSEBISECTION(F, X)   Compute F^{-1}(X) using Bisection.
            % Taken from cumsum as part of chebfun.
            % chebfun/inv.m
            %
            % Copyright 2017 by The University of Oxford and The Chebfun Developers.
            % See http://www.chebfun.org/ for Chebfun information.
            
            a = yMin*ones(size(x));
            b = yMax*ones(size(x));
            c = (a + b)/2;
            
            while ( norm(b - a, inf) >= tol )
                vals = feval(f, c);
                % Bisection:
                I1 = ((vals-x) <= -tol);
                I2 = ((vals-x) >= tol);
                I3 = ~I1 & ~I2;
                a = I1.*c + I2.*a + I3.*c;
                b = I1.*b + I2.*c + I3.*c;
                c = (a+b)/2;
            end
            
            y = c;
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Chebyshev Methods
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [zLobatto, rho_zCheb] = ProjectOntoChebyshevPolynomialsWithTolerance(zIn, rhoFunc, tol)
            m = 3;
            m_max = 15;
            
            zMax = max(zIn);
            zMin = min(zIn);
            Lz = zMax - zMin;
                        
            n = 2^m + 1;
            cutoff = n;
            while (cutoff == n && m <= m_max)
                m = m + 1;
                n = 2^m + 1;
                
                zLobatto = (Lz/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + zMin;
                rho_zCheb = InternalModesSpectral.fct(rhoFunc(zLobatto));
                cutoff = InternalModesSpectral.standardChop(rho_zCheb, tol);
            end
            
            if cutoff < n
                n = cutoff;
                zLobatto = (Lz/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + zMin;
            else
                disp('Unable to project density function within requested tolerance! Using maximum allowed length.');
            end
            
            zLobatto(1) = zMax;
            zLobatto(end) = zMin;
            rho_zCheb = rho_zCheb(1:n);
        end
                
        % Fast Chebyshev Transform
        % By Allan P. Engsig-Karup, apek@imm.dtu.dk.
        function uh = fct(u)
            N  = length(u);
            u  = ifft(u([1:N N-1:-1:2])); % reverse ordering due to Matlab's fft
            uh = ([u(1); 2*u(2:(N-1)); u(N)]);
            
            if any(imag(uh))
                disp('Fast Chebyshev Transform returned imaginary values. Something went wrong!')
            end
            
            uh = real(uh);
        end
        
        % Inverse Fast Chebyshev Transform
        function u = ifct(uh)
            N = length(uh) - 1;
            s = N*[uh(1)*2; uh(2:N); uh(end)*2];
            u = ifft([s; flip(s(2:end-1))],'symmetric');
            u=u(1:N+1);
        end
        
        function bool = IsChebyshevGrid(z_in)
            % make sure the grid is monotonically decreasing
            if (z_in(2) - z_in(1)) > 0
                z_in = flip(z_in);
            end
            
            z_norm = InternalModesSpectral.ChebyshevPolynomialsOnGrid( z_in );
            N_points = length(z_in);
            xi=(0:N_points-1)';
            lobatto_grid = cos(xi*pi/(N_points-1));
            z_diff = z_norm-lobatto_grid;
            if max(abs(z_diff)) < 1e-6
                bool = 1;
            else
                bool = 0;
            end
        end
        
        % Given some Lobatto grid and some desired output grid, return the
        % transformation function T that goes from spectral to the output
        % grid. This basically gives you spectral interpolation.
        function [T, doesOutputGridSpanDomain] = ChebyshevTransformForGrid(lobatto_grid, output_grid)
            if(min(output_grid) < min(lobatto_grid) || max(output_grid) > max(lobatto_grid))
               error('The output grid must be bounded by the lobatto grid'); 
            end
            if (min(output_grid) == min(lobatto_grid) && max(output_grid) == max(lobatto_grid))
                doesOutputGridSpanDomain = 1;
            else
                doesOutputGridSpanDomain = 0;
            end
            
            % T_out transforms vector solutions of the eigenvalue problem
            % into gridded solution on z_out
            if doesOutputGridSpanDomain == 1 && InternalModesSpectral.IsChebyshevGrid(output_grid) == 1
                if length(output_grid) == length(lobatto_grid)
                    T = @(f_cheb) InternalModesSpectral.ifct(f_cheb);
                elseif length(output_grid) > length(lobatto_grid)
                    T = @(f_cheb) InternalModesSpectral.ifct(cat(1,f_cheb,zeros(length(output_grid)-length(lobatto_grid),1)));
                elseif length(output_grid) < length(lobatto_grid)
                    T = @(f_cheb) InternalModesSpectral.ifct(f_cheb(1:length(output_grid)));
                end
            else
                L = max(lobatto_grid)-min(lobatto_grid);
                x = (2/L)*(output_grid-min(lobatto_grid)) - 1;
                t = acos(x);
                TT = zeros(length(output_grid),length(lobatto_grid));
                for iPoly=0:(length(lobatto_grid)-1)
                    TT(:,iPoly+1) = cos(iPoly*t);
                end
                T = @(f_cheb) TT*f_cheb;
            end
        end
        
        function v_p = DifferentiateChebyshevVector(v)
            v_p = zeros(size(v));
            k = length(v)-1;
            v_p(k) = 2*k*v(k+1);
            for k=(length(v)-2):-1:1
                v_p(k) = 2*k*v(k+1) + v_p(k+2);
            end
            v_p(1) = v_p(1)/2;
        end
        
        function s = IntegrateChebyshevVectorWithLimits(v,x,a,b)
            % v are the coefficients of the chebyshev polynomials
            % x is the domain, a Gauss-Lobatto grid
            % a and b are the lower and upper limits, respectively
            v_p = InternalModesSpectral.IntegrateChebyshevVector(v);
            
            s = InternalModesSpectral.ValueOfFunctionAtPointOnGrid(b,x,v_p) - InternalModesSpectral.ValueOfFunctionAtPointOnGrid(a,x,v_p);
        end
                   
        function v_p = IntegrateChebyshevVector(v)
            % Taken from cumsum as part of chebfun.
            % chebfun/@chebtech/cumsum.m
            %
            % Copyright 2017 by The University of Oxford and The Chebfun Developers.
            % See http://www.chebfun.org/ for Chebfun information.
            
            % integration target
            n = length(v);
            v_p = zeros(n+1,1);
            
            % zero-pad coefficients
            v = reshape(v,[],1);
            v = [v; zeros(2,1)];
            
            % Compute b_(2) ... b_(n+1):
            v_p(3:n+1,:) = (v(2:n,:) - v(4:n+2,:)) ./ (2*(2:n).');
            v_p(2,:) = v(1,:) - v(3,:)/2;        % Compute b_1
            t = ones(1, n);
            t(2:2:end) = -1;
            v_p(1,:) = t*v_p(2:end,:);             % Compute b_0 (satisfies f(-1) = 0)
        end
        
        function D = ChebyshevDifferentiationMatrix(n)
            %% Chebyshev Differentiation Matrix
            % Returns the Chebyshev differentiation matrix for the first n polynomials.
            D = zeros(n,n);
            for i=1:n
                for j=1:n
                    if ( j >= i+1 && mod(i+j,2)==1 )
                        D(i,j) = 2*(j-1);
                    else
                        D(i,j) = 0.0;
                    end
                end
            end
            D(1,:)=0.5*D(1,:);
        end
        
        function D = ChebyshevInterpolationDerivative(n)
            %% Chebyshev Interpolation Derivative
            % taken from Canuto, et al. 2.4.33
            D = zeros(n,n);
            N = n-1;
            c = @(j) (j == 0 || j == N)*2 + (j>0 && j<N)*1;
            for j=0:(n-1)
                for l=0:(n-1)
                    if j ~= l
                        D(j+1,l+1) = -(c(j)/c(l))*((-1)^(j+l))/( sin( (j+l)*pi/(2*N) ) * sin( (j-l)*pi/(2*N) ) )/2;
                    elseif j == l && j == 0
                        D(j+1,l+1) = (2*N*N+1)/6;
                    elseif j == l && j == N
                        D(j+1,l+1) = -(2*N*N+1)/6;
                    else
                        D(j+1,l+1) = -cos(pi*j/N)/(2*(sin(j*pi/N))^2);
                    end
                end
            end
            D(1,:)=D(1,:);
        end
        
        function value = ValueOfFunctionAtPointOnGrid( x0, x, func_cheb )
           % We have the Chebyshev coefficents of function func_cheb, defined on grid x, return the value at x0;
           L = max(x)-min(x);
           x_norm = (2/L)*(x0-min(x)) - 1;
           t = acos(x_norm);
           
           N_polys = length(func_cheb);
%            value=sum(func_cheb.*cos(t*(0:(N_polys-1))));   
           value = func_cheb(1)*ones(size(x0));
           for i=2:N_polys
               value = value + func_cheb(i)*cos(t*(i-1));
           end
                   
        end
        
        function [varargout] = ChebyshevPolynomialsOnGrid( x, N_polys )
            %% Chebyshev Polynomials on Grid
            % Compute the the first N Chebyshev polynomials and their derivatives for
            % an arbitrary grid x.
            %
            % x_norm = ChebyshevPolynomialsOnGrid( x ) with exactly one argument, x,
            % returns the x normalized to its typical [-1,1] values.
            %
            % T = ChebyshevPolynomialsOnGrid( x, N_polys ) returns the first N_poly
            % Chebyshev polynomials for an arbitrary grid x.
            %
            % [T, T_x, T_xx,...] = ChebyshevPolynomialsOnGrid( x, N_polys ) returns the
            % first N_poly Chebyshev polynomials and their derivatives for an arbitrary
            % grid x.
            %
            % The returned matrices T, T_xx, etc are size(T) = [length(x) N_polys],
            % i.e., the polynomials are given column-wise.
            
            N_points = length(x);
            
            % These are the normalized coordinates for Chebyshev polynomials.
            L = max(x)-min(x);
            x_norm = (2/L)*(x-min(x)) - 1;
            t = acos(x_norm);
            
            % if there's only one input argument, we just return x_norm
            if nargin == 1
                varargout{1} = x_norm;
                return;
            else
                if N_polys < 4
                    disp('You must request at least four polynomials. Fixing that for you.')
                    N_polys = 4;
                end
            end
            
            N_diff = nargout-1;
            varargout = cell(1,nargout);
            
            % It's easy to create the polynomials, they're stretched cosines!
            T = zeros(N_points,N_polys);
            for iPoly=0:(N_polys-1)
                T(:,iPoly+1) = cos(iPoly*t);
            end
            
            varargout{1} = T;
            
            % Now use the recursion formula to compute derivates of polynomials.
            for n=1:N_diff
                T = varargout{n};
                T_x = zeros(size(T));
                T_x(:,2) = T(:,1);
                T_x(:,3) = 2*2*T(:,2);
                for j=4:N_polys
                    m = j-1;
                    T_x(:,j) = (m/(m-2))*T_x(:,j-2) + 2*m*T(:,j-1);
                end
                varargout{n+1} = (2/L)*T_x;
            end
            
        end
        
        function [z_lobatto_grid] = FindSmallestChebyshevGridWithNoGaps(z)
            % Want to create a chebyshev grid that never has two or more point between
            % its points. If that makes sense.
            if (z(2) - z(1)) > 0 % make z_out decreasing
                z = flip(z);
            end
            
            L = max(z)-min(z);
            np = ceil(log2(length(z)));
            n = 2^np;
            z_lobatto_grid = (L/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + min(z);
            
            while( length(unique(interp1(z_lobatto_grid,z_lobatto_grid,z,'previous'))) ~= length(z) )
                np = np + 1;
                n = 2^np;
                z_lobatto_grid = (L/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + min(z);
            end
        end
        
        function cutoff = standardChop(coeffs, tol)
            % Copyright 2017 by The University of Oxford and The Chebfun Developers. 
            % See http://www.chebfun.org/ for Chebfun information.
            %
            % This is taken, without comments and safe checks, from the
            % above developers. They get sole credit. Jared Aurentz and Nick Trefethen, July 2015.
            n = length(coeffs);
            cutoff = n;
            if ( n < 17 )
                return
            end
            
            envelope = cummax(abs(coeffs),'reverse');
            if envelope(1) == 0
                cutoff = 1;
                return
            else
                envelope = envelope/envelope(1);
            end
            
            for j = 2:n
                j2 = round(1.25*j + 5);
                if ( j2 > n )
                    % there is no plateau: exit
                    return
                end
                e1 = envelope(j);
                e2 = envelope(j2);
                r = 3*(1 - log(e1)/log(tol));
                plateau = (e1 == 0) | (e2/e1 > r);
                if ( plateau )
                    % a plateau has been found: go to Step 3
                    plateauPoint = j - 1;
                    break
                end
            end
            
            
            if ( envelope(plateauPoint) == 0 )
                cutoff = plateauPoint;
            else
                j3 = sum(envelope >= tol^(7/6));
                if ( j3 < j2 )
                    j2 = j3 + 1;
                    envelope(j2) = tol^(7/6);
                end
                cc = log10(envelope(1:j2));
                cc = cc(:);
                cc = cc + linspace(0, (-1/3)*log10(tol), j2)';
                [~, d] = min(cc);
                cutoff = max(d - 1, 1);
            end
            
        end
        
        function roots = FindRootsFromChebyshevVector(f_cheb, zLobatto)
            % Copyright (c) 2007, Stephen Morris 
            % All rights reserved.
            n = length(f_cheb);
            
            f_cheb(abs(f_cheb)<1e-12) = 1e-15;
            
            A=zeros(n-1);   % "n-1" because Boyd's notation includes the zero-indexed
            A(1,2)=1;       % elements whereas Matlab's of course does not allow this.
            % In other words, to replicate Boyd's N=5 matrix we need to
            % set n=6.
            for j=2:n-2
                for k=1:n-1
                    if j==k+1 || j==k-1
                        A(j,k)=0.5;
                    end
                end
            end
            for k=1:n-1
                A(n-1,k)=-f_cheb(k)/(2*f_cheb(n));  % c(1) in our notation is c(0) in Boyd's
            end
            A(n-1,n-2)=A(n-1,n-2)+0.5;
            % Now we have the companion matrix, we can find its eigenvalues using the
            % MATLAB built-in function.
            eigvals=eig(A);
            
            % We're only interested in the real elements of the matrix:
            realvals=(arrayfun(@(x) ~any(imag(x)),eigvals)).*eigvals;
            
            % Of course these are the roots scaled to the canonical interval [-1,1]. We
            % need to map them back onto the interval [a,b]; we widen the interval just
            % a fraction to make sure that we don't miss any that are right on the
            % edge.
            a = min(zLobatto);
            b = max(zLobatto);
            rangevals=nonzeros((arrayfun(@(x) abs(x)<=1.001, realvals)).*realvals);
            roots=sort((rangevals.*0.5*(b-a)) + (0.5*(b+a)));
        end
        
    end
end


