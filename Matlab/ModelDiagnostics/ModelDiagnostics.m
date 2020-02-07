classdef ModelDiagnostics < handle
    % ModelDiagnostics
    properties (Access = public)
        Lx, Ly, Lz % Domain size
        Nx, Ny, Nz % Number of points in each direction
        latitude
        
        x, y, z
        k, l, j
        X,Y,Z
        K,L,J
        
        rhoBar % mean density, function of z
        rho0 % density at the surface
        
        N2 % will be forced to be 1x1xLz
        Nmax % sqrt(max(N2))
        K2, Kh, f0
        
        u,v,w,rhoPrime
    end
    
    properties (Dependent)
        % Alternative forms of density background and density anomaly
        bbar        % bbar = -(g/rho_0)*rhobar
        b           % buoyancy anomaly, b = -(g/rho_0)*rhoPrime
        eta         % eta = - rhoPrime/rhobar_z
        
        
        bbar_z      % bbar_z = N^2 = -(g/rho_0)*rhobar_z
        rhoBar_z    % rhobar_z = -(rho0/g)*N^2
        
        % Single derivatives of the primary variables
        u_x, u_y, u_z
        v_x, v_y, v_z
        w_x, w_y, w_z
        b_x, b_y, b_z
        eta_x, eta_y, eta_z
        
        % Components of vorticity
        vorticity_x	% w_y - v_z
        vorticity_y	% u_z - w_x
        vorticity_z	% v_x - u_y
    end
    
    properties (Constant)
        g = 9.81;
        L_gm = 1.3e3; % thermocline exponential scale, meters
        N_gm = 5.2e-3; % reference buoyancy frequency, radians/seconds
        E_gm = 6.3e-5; % non-dimensional energy parameter
        E = (1.3e3)*(1.3e3)*(1.3e3)*(5.2e-3)*(5.2e-3)*(6.3e-5);
    end
    
    methods
        function self = ModelDiagnostics(dims, n, latitude, rhoBar, N2)
            if length(dims) ~=3 || length(n) ~= 3
                error('The dims and n variables must be of length 3. You need to specify x,y,z');
            end
            
            if length(rhoBar) ~= dims(3) || length(N2) ~= dims(3)
               error('N2 must have the same number of points as the z dimension'); 
            end
            
            self.Lx = dims(1);
            self.Ly = dims(2);
            self.Lz = dims(3);
            
            self.Nx = n(1);
            self.Ny = n(2);
            self.Nz = n(3);
            
            self.latitude = latitude;
            
            dx = self.Lx/self.Nx;
            dy = self.Ly/self.Ny;
            self.x = dx*(0:self.Nx-1)'; % periodic basis
            self.y = dy*(0:self.Ny-1)'; % periodic basis
            self.z = z; % cosine basis (not your usual dct basis, however)
            
            self.rhoBar = rhoBar;
            self.rho0 = rhoBar(end);
            
            % Let's make N2 exist in the 3rd dimension, for east
            % multiplication.
            if ~exist('N2','var')
                N2 = -(self.g/self.rho0)*DiffCosine(self.z,rhoBar);
            end
            self.N2 = reshape(N2,1,1,[]);
            self.Nmax = sqrt(max(self.N2));
            
            % Spectral domain, in radians
            if self.Nx > 0 && self.Lx > 0
                dk = 1/self.Lx;          % fourier frequency
                self.k = 2*pi*([0:ceil(self.Nx/2)-1 -floor(self.Nx/2):-1]*dk)';
            else
                self.k = 0;
            end
            
            if self.Ny > 0 && self.Ly > 0
                dl = 1/self.Ly;          % fourier frequency
                self.l = 2*pi*([0:ceil(self.Ny/2)-1 -floor(self.Ny/2):-1]*dl)';
            else
                self.l = 0;
            end
            
            nModes = length(z);
            self.j = (1:nModes)';
            
            [K,L,J] = ndgrid(self.k,self.l,self.j);
            [X,Y,Z] = ndgrid(self.x,self.y,self.z);
            
            self.L = L; self.K = K; self.J = J;
            self.X = X; self.Y = Y; self.Z = Z;
            
            self.f0 = 2 * 7.2921E-5 * sin( self.latitude*pi/180 );
            self.K2 = self.K.*self.K + self.L.*self.L;   % Square of the horizontal wavenumber
            self.Kh = sqrt(self.K2);
        end
         
        function InitializeWithHorizontalVelocityAndDensityPerturbationFields(self, u, v, w, rhoPrime)
            self.u = u;
            self.v = v;
            self.w = w;
            self.rhoPrime = rhoPrime;
        end
        
        function InitializeWithHorizontalVelocityAndIsopycnalDisplacementFields(self, u, v, w, eta)
            self.u = u;
            self.v = v;
            self.w = w;
            self.rhoPrime = -eta .* self.rhoBar_z;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %  Other ways of expressing the density perturbation
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function value = get.bbar(self)
            value = -(self.g/self.rho0)*self.rhoBar;
        end
        
        function value = get.b(self)
            value = -(self.g/self.rho0)*self.rhoPrime;
        end
        
        function value = get.eta(self)
            value = -self.rhoPrime./self.rhoBar_z;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %  Other ways of expressing the background density gradient
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function value = get.bbar_z(self)
            value = self.N2;
        end
        
        function value = get.rhoBar_z(self)
            value = -(self.rho0/self.g)*self.N2;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %  Single derivatives of the primary variables
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function value = get.u_x(self)
            value = DiffFourier(self.x,self.u,1,1);
        end
        
        function value = get.u_y(self)
            value = DiffFourier(self.y,self.u,1,2);
        end
        
        function value = get.u_z(self)
            value = DiffCosine(self.z,self.u);
        end
        
        function value = get.v_x(self)
            value = DiffFourier(self.x,self.v,1,1);
        end
        
        function value = get.v_y(self)
            value = DiffFourier(self.y,self.v,1,2);
        end
        
        function value = get.v_z(self)
            value = DiffCosine(self.z,self.v);
        end
        
        function value = get.w_x(self)
            value = DiffFourier(self.x,self.w,1,1);
        end
        
        function value = get.w_y(self)
            value = DiffFourier(self.y,self.w,1,2);
        end
        
        function value = get.w_z(self)
            value = DiffSine(self.z,self.w);
        end
        
        function value = get.b_x(self)
            value = DiffFourier(self.x,self.b,1,1);
        end
        
        function value = get.b_y(self)
            value = DiffFourier(self.y,self.b,1,2);
        end
        
        function value = get.b_z(self)
            value = DiffSine(self.z,self.b);
        end
        
        function value = get.eta_x(self)
            value = DiffFourier(self.x,self.eta,1,1);
        end
        
        function value = get.eta_y(self)
            value = DiffFourier(self.y,self.eta,1,2);
        end
        
        function value = get.eta_z(self)
            value = DiffSine(self.z,self.eta);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %  Useful derived quantities
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function value = get.vorticity_x(self)
            value = self.w_y - self.v_z; % w_y - v_z
        end
        
        function value = get.vorticity_y(self)
            value = self.u_z - self.w_x; % u_z - w_x
        end
        
        function value = get.vorticity_z(self)
            value = self.v_x - self.u_y; % v_x - u_y
        end
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Measures of Potential Vorticity
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        function PV = ErtelPV(self,noBackground)
            %% Ertel PV,  1/s^3
            % PV = (N^2(z) \hat{k} + \nabla b) \cdot (f \hat{k} + \omega)
            %
            % The background PV, f0*N2, can be safely eliminated in constant
            % stratification, but it remains important in variable
            % stratification due to vertical advection.
            
            PVx = self.vorticity_x .* self.b_x;
            PVy = self.vorticity_y .* self.b_y;
            if exist('noBackground','var') && noBackground == 1
                PVz = (self.f0 + self.vorticity_z) .* self.b_z + self.vorticity_z .* self.bbar_z;
            else
                PVz = (self.f0 + self.vorticity_z) .* (self.b_z + self.bbar_z);
            end
            
            PV = PVx + PVy + PVz;
        end
        
        function PV = LinearErtelPV(self,noBackground)
            %% Linear Ertel PV  (the linear terms in Ertel PV) 1/s^3
            % PV = f_0 b_z + (f0 + vorticity_z) *N^2(z)
            %
            % Same as above, the background PV, f0*N2, can be safely
            % eliminated in constant stratification, but it remains
            % important in variable stratification due to vertical
            % advection.
            
            if exist('noBackground','var') && noBackground == 1
                PV = self.f0 * self.b_z + self.vorticity_z .* self.N2;
            else
                PV = self.f0 * self.b_z + (self.f0 + self.vorticity_z) .* self.N2;
            end
        end
        
        function PV = QGPV(self)
            %% QG PV, quasigeostrophic potential vorticity
            % PV = N^2 (-f_0 eta_z + \zeta_z + f0) where eta = -rho/rhobar_z
            %
            % This is only the same as linear Ertel PV when N^2 is
            % constant.
            %
            % We include the f0 term to keep it comparable in magnitude to
            % Ertel PV.
            
            PV = self.N2 .* ( self.vorticity_z - self.f0 * self.eta_z );
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Variances and energy as a function of depth
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [Euv,z] = HorizontalVelocityVariance(self)
            %% Horizontal velocity variance as a function of depth, m^2/s^2
            %
            Euv = squeeze(mean(mean((self.u).^2 + (self.v).^2,1),2));
            z = self.z;
        end
        
        function [Ew,z] = VerticalVelocityVariance(self)
            %% Vertical velocity variance as a function of depth, m^2/s^2
            %
            Ew = squeeze(mean(mean((self.w).^2,1),2));
            z = self.z;
        end
        
        function [Eeta,z] = IsopycnalVariance(self)
            %% Isopycnal velocity variance as a function of depth, m^2
            %
            Eeta = squeeze(mean(mean((self.eta).^2,1),2));
            z = self.z;
        end
        
        function [Euv,z] = HorizontalVelocityVarianceWKB(self)
            %% WKB scaled horizontal velocity variance as a function of depth, m^2/s^2
            %
            Euv = (self.N_gm./sqrt(self.N2)) .* self.HorizontalVelocityVariance();
            z = self.z;
        end
        
        function [Ew,z] = VerticalVelocityVarianceWKB(self)
            %% WKB scaled vertical velocity variance as a function of depth, m^2/s^2
            %
            Ew = (sqrt(self.N2)./self.N_gm) .* self.VerticalVelocityVariance();
            z = self.z;
        end
        
        function [Eeta,z] = IsopycnalVarianceWKB(self)
            %% WKB scaled sopycnal velocity variance as a function of depth, m^2
            %
            Eeta = (sqrt(self.N2)./self.N_gm) .* self.IsopycnalVariance();
            z = self.z;
        end
        
    end
end