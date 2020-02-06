classdef ModelDiagnostics < handle
    % ModelDiagnostics
    properties (Access = public)
        Lx, Ly, Lz % Domain size
        Nx, Ny, Nz % Number of points in each direction
        nModes
        latitude
        
        x, y, z
        k, l, j
        X,Y,Z
        K,L,J
        
        N2
        K2, Kh, f0
    end
    
    methods
        function self = ModelDiagnostics(dims, n, N2, latitude)
            if length(dims) ~=3 || length(n) ~= 3
                error('The dims and n variables must be of length 3. You need to specify x,y,z');
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
            
            self.N2 = N2;
            
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
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % 
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        function PV = ErtelPV(self, u, v, w, rho_prime)
            %% Ertel PV
            % PV = (N^2(z) \hat{k} + \nabla b) \cdot (f \hat{k} + \omega)
        end
        
        function PV = LinearErtelPV(self, u, v, w, rho_prime)
            %% Linear Ertel PV ? linear terms in Ertel PV
            % PV = f_0 b_z + \zeta_z N^2(z)
        end
        
        function PV = QGPV(self, u, v, w, rho_prime)
            %% QG PV ? Quasigeostrophic potential vorticity
            % This is only the same as linear Ertel PV when N^2 is
            % constant.
            % PV = N^2 (-f_0 eta_z + \zeta_z ) where eta = -rho/rhobar_z
            
        end
        
        function [ErtelPV, LinearErtelPV, QGPV] = PV(self, u, v, w, rho_prime)
            %% Ertel PV
            % PV = (N^2(z) \hat{k} + \nabla b) \cdot (f \hat{k} + \omega)
            zeta_x = DiffFourier(self.y,w,1,2) - DiffCosine(self.z,v); % w_y - v_z
            zeta_y = DiffCosine(self.z,u) - DiffFourier(self.x,w,1,1); % u_z - w_x
            zeta_z = DiffFourier(self.x,v,1,1) - DiffFourier(self.y,u,1,2); % v_x - u_y
            
            % scaling b so that it is meters (isopycnal height), with a sign difference
            b = -(wavemodel.g/wavemodel.rho0)*rho_prime/N0/N0;
            
            % The derivatives are unitless
            b_x = DiffFourier(self.x,b,1,1);
            b_y = DiffFourier(self.y,b,1,2);
            b_z = DiffSine(self.z,b);
            
            % this should be 1, if we've done this correctly
            bbar_z = wavemodel.N2AtDepth(wavemodel.Z)/N0/N0;
            
            PV_x = zeta_x .* b_x;
            PV_y = zeta_y .* b_y;
            PV_z = (zeta_z + self.f0) .* b_z + zeta_z .* bbar_z;
            
            LinearErtelPV = zeta_z .* bbar_z + self.f0 * b_z;
            ErtelPV = PV_x + PV_y + PV_z;
            
        end
        
    end
end