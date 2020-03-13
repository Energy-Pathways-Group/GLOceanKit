classdef Boussinesq2D < handle
    % Boussinesq2D Nonlinear model
    %
    % Jeffrey J. Early
    % jeffrey@jeffreyearly.com
    %
    % March 12th, 2020      Version 1.0

    properties (Access = public)
        Lx, Lz % Domain size
        Nx, Nz % Number of points in each direction
        
        x, z
        k, m, m_s
        X,Z
        K,M,M_s
        
        rhobar, rho0
        N2, Nmax
        
        version = 1.0
        
        nu_x, nu_z
        
        psi_n, b_n, nabla2_psi_n
        
        shouldAntialias = 0;
    end
        
    properties (Constant)
        g = 9.81;
    end
    
%     methods(Abstract, Access = public)
%         N2 = N2AtDepth(self,z)
%         rho = RhoBarAtDepth(self,z)
%     end

    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = Boussinesq2D(dims, n, z, N2)
            if length(dims) ~=2 || length(n) ~= 2
                error('The dims and n variables must be of length 2. You need to specify x,z');
            end
            
            self.Lx = dims(1);
            self.Lz = dims(2);
            
            self.Nx = n(1);
            self.Nz = n(2);
                       
            dx = self.Lx/self.Nx;            
            self.x = dx*(0:self.Nx-1)'; % periodic basis
            self.z = z; % cosine dct-I basis
            
            self.N2 = N2;
            self.Nmax = sqrt(max(N2));
%             self.rhobar = self.RhoBarAtDepth(self.z);
            
            % Using a Fast Fourier Transform
            dk = 1/self.Lx;          % fourier frequency
            self.k = 2*pi*([0:ceil(self.Nx/2)-1 -floor(self.Nx/2):-1]*dk)';
            
            % Using a Fast Discrete Cosine Transform (DCT-I)
            dm = 2*pi/(2*(self.Nz-1)*(z(2)-z(1)));
            self.m = dm*(0:(self.Nz-1))';
            self.m_s = self.m(2:end-1);
            
            [self.K,self.M] = ndgrid(self.k,self.m);
            [self.K,self.M_s] = ndgrid(self.k,self.m_s);
            [self.X,self.Z] = ndgrid(self.x,self.z);                                    
        end
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Create a single wave (public)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [omega,h] = InitializeWithPlaneWave(self, k0, j0, U)
            im = InternalModesConstantStratification(self.Nmax,[min(self.z) max(self.z)],self.z,0);
            im.normalization = Normalization.uMax;
            [~,G,h,omega] = im.ModesAtWavenumber(k0);
            
            self.psi_n = U*h(j0)*cos(k0*self.X).*(G(:,j0).');
            self.b_n = -(U*k0*h(j0)/omega(j0))*cos(k0*self.X).*(G(:,j0).');
            h = h(j0);
            omega = omega(j0);
            
            self.nabla2_psi_n = self.nabla2_psi(self.psi_n);
        end
        
        function [f_psi, f_b] = linearFlux(nabla2_psi,b)
            nabla2_psi_bar = TransformForwardFS(nabla2_psi);
            b_bar = TransformForwardFS(b);
            
            f_psi = -self.b_x(b_bar);
            f_b = -self.N2*self.psi_x(nabla2_psi_bar);
        end
        
        function f = linearFluxCat(y0)
            nabla2_psi = y0(:,:,1);
            b = y0(:,:,2);
            
            nabla2_psi_bar = TransformForwardFS(nabla2_psi);
            b_bar = TransformForwardFS(b);
            
            f_psi = -self.b_x(b_bar);
            f_b = -self.N2*self.psi_x(nabla2_psi_bar);
            
            f = cat(3,f_psi,f_b);
        end
        
        function Q = SVV(self)
            % Builds the spectral vanishing viscosity operator
            k_max = max(self.k);
            m_max = max(self.m);
            if self.shouldAntialias == 1
                k_max = 2*k_max/3;
                m_max = 2*m_max/3;
            end
            
            dk = self.k(2)-self.k(1);
            dm = self.m(2)-self.m(1);
            k_cutoff = dk*(k_max/dk)^(3/4);
            m_cutoff = dm*(m_max/dm)^(3/4);
            
            Q = exp( - ((self.K-k_max)./(self.K-k_cutoff)).^2 );
            Q = Q.*exp( - ((self.M-m_max)./(self.M-m_cutoff)).^2 );
            Q(self.K<k_cutoff & self.M<m_cutoff) = 0;
            Q(self.K>k_max & self.M>m_cutoff) = 1;
        end
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Transformations between spatial and spectral domains
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function psi_bar = TransformForwardFS(self,psi)
            % Transform from physical coordinates to a (fourier,sine) basis
            psi_bar = SineTransformForward(self.z,psi,2,'both');
            psi_bar = FourierTransformForward(self.x,psi_bar,1);
        end
        
        function psi = TransformBackFS(self,psi_bar)
            % Transform from (fourier,sine) basis to physical coordinates
            psi_bar = FourierTransformBack( self.k, psi_bar, 1 );
            psi = SineTransformBack( self.m_s, psi_bar, 2 );
        end
        
        function psi = TransformBackFC(self,psi_bar)
            % Transform from (fourier,cosine) basis to physical coordinates
            psi_bar = cat(2,zeros(self.Nx,1),psi_bar,zeros(self.Nx,1));
            psi_bar = FourierTransformBack( self.k, psi_bar, 1 );
            psi = CosineTransformBack( self.m, psi_bar, 2); 
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Differential operators
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function damp_psi = damp_psi(self,nabla2_psi_bar)
            L = - self.SVV .* (self.nu_x*self.K.^4 + self.nu_z*self.M_s.^4)./(self.K.* self.K + self.M_s.*self.M_s);
            damp_psi = self.TransformBackFS(L.*nabla2_psi_bar);
        end
        
        function psi_bar = AntiAlias(self,psi_bar)
            % Should anti-alias *after* each time step, before transforming
            % back into spectral space (and before advecting particles!)
            AA = ones(size(self.K));
            AA( sqrt(self.K.*self.K/max(self.k)^2 + self.M_s.*self.M_s/max(self.m)^2) > 2/3 ) = 0;
            psi_bar = AA .* psi_bar;
        end
        
        function nabla2_psi = nabla2_psi(self,psi)
            psi_bar = self.TransformForwardFS( psi );
            nabla2_psi = self.TransformBackFS( -(self.K.* self.K + self.M_s.*self.M_s).*psi_bar );
        end
        
        function psi_x = psi_x(self,nabla2_psi_bar)
            L = sqrt(-1)*self.K./(self.K.* self.K + self.M_s.*self.M_s);
            psi_x = self.TransformBackFS( L .* nabla2_psi_bar );
        end
        
        function psi_z = psi_z(self,nabla2_psi_bar)
            L = self.M_s./(self.K.* self.K + self.M_s.*self.M_s);
            psi_z = self.TransformBackFC( L .* nabla2_psi_bar );
        end
        
        function nabla2_psi_x = nabla2_psi_x(self,nabla2_psi_bar)
            L = sqrt(-1)*self.K;
            nabla2_psi_x = self.TransformBackFS( L .* nabla2_psi_bar );
        end
        
        function nabla2_psi_z = nabla2_psi_z(self,nabla2_psi_bar)
            L = self.M_s;
            nabla2_psi_z = self.TransformBackFC( L .* nabla2_psi_bar );
        end
        
        function b_x = b_x(self,b)
            b_x = DiffFourier(self.x,b,1,1);
        end
        
        function b_z = b_z(self,b)
            b_z = DiffSine(self.z,b,1,2);
        end
        
    end
    
    % For the psi equation, we need:
    % (dx) psi
    % (dz) \nabla^2 psi
    % (dz) psi
    % (dx) \nabla^2 psi
    % (dx) b
    %
    % For the b equation, we need:
    % (dx) psi
    % (dz) b
    % (dz) psi
    % (dx) b
    %
    % Operators needed
    %
    % \nabla^2 psi -> psi_x, so dx*\nabla^{-2}
    % \nabla^2 psi -> psi_z, so dz*\nabla^{-2}
    % \nabla^2 psi -> nabla^2 psi_z, so dz
    % \nabla^2 psi -> nabla^2 psi_x, so dx
    % 
    
    methods (Static)
        
    end
    
end


