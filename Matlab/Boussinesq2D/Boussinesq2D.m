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
        
        x,z
        k,m,m_s
        X,Z
        K,M_s
        
        rhobar, rho0
        N2
        internalModes
        
        version = 1.0
        
        nu_x, nu_z
                
        dt
        integrator
        t = 0
        y % cell array with {nabla2_psi,b,xi,zeta}
        
        nParticles = 0
        shouldAntialias = 0;
        nonlinear = 1; % set to 0 to evolve the equations linearly
    end
        
    properties (Constant)
        g = 9.81;
    end
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = Boussinesq2D(dims, n, z, rhobar)
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
            
            self.internalModes = InternalModes(rhobar,[min(self.z) max(self.z)],self.z,0);
            
            self.N2 = self.internalModes.N2.';
            
            % Using a Fast Fourier Transform
            dk = 1/self.Lx;          % fourier frequency
            self.k = 2*pi*([0:ceil(self.Nx/2)-1 -floor(self.Nx/2):-1]*dk)';
            
            % Using a Fast Discrete Cosine Transform (DCT-I)
            dm = 2*pi/(2*(self.Nz-1)*(z(2)-z(1)));
            self.m = dm*(0:(self.Nz-1))';
            self.m_s = self.m(2:end-1);
                        
            [self.K,self.M_s] = ndgrid(self.k,self.m_s);
            [self.X,self.Z] = ndgrid(self.x,self.z);
            
            self.y = {zeros(size(self.K)), zeros(size(self.K)),[],[]};
        end
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Create a single wave (public)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [omega,h] = InitializeWithPlaneWaveCorrection(self, k0, j0, U)
            self.internalModes.normalization = Normalization.uMax;
            [~,G,h,omega] = self.internalModes.ModesAtWavenumber(k0);
            
            h = h(j0);
            omega = omega(j0);
            G = G(:,j0).';
            c2 = self.g * h;
            
            psi_n = U*h*cos(k0*self.X).*G;
            b_n = self.N2.*((U*k0*h/omega)*cos(k0*self.X)).*G;
            
            self.internalModes.normalization = Normalization.uMax;
            [~,G2k,h2k] = self.internalModes.ModesAtWavenumber(2*k0);
%             N = InternalModes.NumberOfWellConditionedModes(G2k);
            ratio = self.internalModes.internalModes.rho_zz ./ self.internalModes.internalModes.rho_z;
            
            f = ratio.*(G.*G).';
            coeffs = G2k\f;
            Gamma = (G2k * ( (h*h2k./(h-h2k)).' .* coeffs)).';
            Gamma_b = 0.5*f.' + Gamma;
            
%             figure
%             plot([G;Gamma;Gamma_b],self.z)
            
            psi1_n = (U*U*h/(2*sqrt(c2)))*cos(2*k0*self.X).*Gamma;
            b1_n = (U*U*h/(2*c2))*self.N2.*cos(2*k0*self.X).*Gamma_b;
            
            self.InitializeWithPsiAndB(psi_n+psi1_n,b_n+b1_n,omega);
        end
        
        function [omega,h] = InitializeWithPlaneWaveNewCorrection(self, k0, j0, U)
            self.internalModes.normalization = Normalization.uMax;
            [~,G,h,omega] = self.internalModes.ModesAtWavenumber(k0);
            
            h = h(j0);
            omega = omega(j0);
            G = G(:,j0).';
            c2 = self.g * h;
            
            psi_n = U*h*cos(k0*self.X).*G;
            b_n = self.N2.*((U*k0*h/omega)*cos(k0*self.X)).*G;
            
            self.internalModes.normalization = Normalization.uMax;
            [~,G2k,h2k] = self.internalModes.ModesAtWavenumber(2*k0);
            %             N = InternalModes.NumberOfWellConditionedModes(G2k);
            ratio = self.internalModes.internalModes.rho_zz ./ self.internalModes.internalModes.rho_z;
            
            f = ratio.*(G.*G).';
            coeffs = G2k\f;
            Gamma = (G2k * ( (h2k./(h-h2k)).' .* coeffs)).';
            Gamma_b = (G2k * ( ((h+h2k)./(h-h2k)).' .* coeffs)).';
            
            %             figure
            %             plot([G;Gamma;Gamma_b],self.z)
            
            psi1_n = (U*U*h*h/(2*sqrt(c2)))*cos(2*k0*self.X).*Gamma;
            b1_n = (U*U*h*h/(4*c2))*self.N2.*cos(2*k0*self.X).*Gamma_b;
            
            self.InitializeWithPsiAndB(psi_n+psi1_n,b_n+b1_n,omega);
        end
        
        function [omega,h] = InitializeWithPlaneWave(self, k0, j0, U)
            self.internalModes.normalization = Normalization.uMax;
            [~,G,h,omega] = self.internalModes.ModesAtWavenumber(k0);
            
            psi_n = U*h(j0)*cos(k0*self.X).*(G(:,j0).');
            b_n = self.N2.*((U*k0*h(j0)/omega(j0))*cos(k0*self.X)).*(G(:,j0).');
            h = h(j0);
            omega = omega(j0);
            
            self.InitializeWithPsiAndB(psi_n,b_n,omega);
        end
        
        function InitializeWithPsiAndB(self,psi_n,b_n,omega)
            % other stuff that needs to be initialized...
            nabla2_psi_n = self.Nabla2PsiFromPsi(psi_n);
            self.y = {nabla2_psi_n;b_n;[];[]};
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Set the viscosity
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            u = self.u;
            w = self.w;
            U = max(u(:));
            W = max(w(:));
            self.nu_x = U*(self.x(2)-self.x(1));
            self.nu_z = W*(self.z(2)-self.z(1));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Set the time-step
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cfl = 0.25;
            dt_u = cfl*(self.x(2)-self.x(1))/U;
            dt_w = cfl*(self.z(2)-self.z(1))/W;
            fprintf('cfl condition for (u,w)=>dt=(%.1f,%.1f) seconds for period of %.1f seconds\n',dt_u,dt_w,2*pi/omega);

            self.dt = min(dt_u,dt_w);
            
            if ~isempty(omega)
                n = round(2*pi/omega/self.dt);
                self.dt = 2*pi/omega/n;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Set up the integrator
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            f = @(t,y0) self.fluxWithParticles(y0);
            
            self.integrator = ArrayIntegrator(f,self.y,self.dt);
        end
        
        function StepForwardToTime(self,time)
            self.y = self.integrator.StepForwardToTime(time);
            self.t = self.integrator.currentTime;
        end
        
        function setParticlePositions(self,xi0,zeta0)
            self.nParticles = length(xi0);
            self.y{3} = xi0;
            self.y{4} = zeta0;
            self.integrator = ArrayIntegrator(@(t,y0) self.fluxWithParticles(y0),self.y,self.dt);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Accessors
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function u = u(self)
            nabla2_psi_bar = self.TransformForwardFS(self.y{1});
            u = self.Psi_zFromNabla2PsiBar(nabla2_psi_bar);
        end
        
        function w = w(self)
            nabla2_psi_bar = self.TransformForwardFS(self.y{1});
            w = -self.Psi_xFromNabla2PsiBar(nabla2_psi_bar);
        end
        
        function [u,w] = uw(self)
            nabla2_psi_bar = self.TransformForwardFS(self.y{1});
            u = self.Psi_zFromNabla2PsiBar(nabla2_psi_bar);
            w = -self.Psi_xFromNabla2PsiBar(nabla2_psi_bar);
        end
        
        function psi = psi(self)
           psi = self.PsiFromNabla2Psi(self.y{1}); 
        end
        
        function b = b(self)
            b = self.y{2};
        end
        
        function xi = xi(self)
            xi = self.y{3};
        end
        
        function zeta = zeta(self)
            zeta = self.y{4};
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % The flux (F) in the equation dy/dt = F(t,y)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
        function f = fluxWithParticles(self,y0)
            f = cell(4,1);
            nabla2_psi = y0{1};
            b = y0{2};
            xi = y0{3};
            zeta = y0{4};

            if self.nonlinear == 1
                nabla2_psi_bar = self.AntiAlias( self.TransformForwardFS(nabla2_psi) );
                u = self.Psi_zFromNabla2PsiBar(nabla2_psi_bar);
                w = -self.Psi_xFromNabla2PsiBar(nabla2_psi_bar);
                b_x = self.b_x(b);
                b_z = self.b_z(b);
                nabla2_psi_x = self.Nabla2Psi_xFromNabla2PsiBar(nabla2_psi_bar);
                nabla2_psi_z = self.Nabla2Psi_zFromNabla2PsiBar(nabla2_psi_bar);
                
                f{1} = -u.*nabla2_psi_x - w.*nabla2_psi_z - b_x; % + self.damp_psi(nabla2_psi_bar);
                f{2} =-u.*b_x - w.*(self.N2 + b_z);
            else
                nabla2_psi_bar = self.TransformForwardFS(nabla2_psi);
                w = -self.Psi_xFromNabla2PsiBar(nabla2_psi_bar);
                
                f{1} = -self.b_x(b); % + self.damp_psi(nabla2_psi_bar);
                f{2} = -self.N2.*w;
            end
            
            if self.nParticles > 0    
                f{3} = interpn(self.X,self.Z,u,xi,zeta);
                f{4} = interpn(self.X,self.Z,w,xi,zeta);
            else
                f{3} = [];
                f{4} = [];
            end
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
%             psi = cat(2,zeros(self.Nx,1),psi,zeros(self.Nx,1));
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
        
        function psi_bar = AntiAlias(self,psi_bar)
            % Should anti-alias *after* each time step, before transforming
            % back into spectral space (and before advecting particles!)
            AA = ones(size(self.K));
            AA( sqrt(self.K.*self.K/max(self.k)^2 + self.M_s.*self.M_s/max(self.m)^2) > 2/3 ) = 0;
            psi_bar = AA .* psi_bar;
        end
        
        function psi = PsiFromNabla2Psi(self,nabla2_psi)
            L = -1./(self.K.* self.K + self.M_s.*self.M_s);
            nabla2_psi_bar = self.TransformForwardFS( nabla2_psi );
            psi = self.TransformBackFS( L .* nabla2_psi_bar );
        end
        
        function nabla2_psi = Nabla2PsiFromPsi(self,psi)
            psi_bar = self.TransformForwardFS( psi );
            nabla2_psi = self.TransformBackFS( -(self.K.* self.K + self.M_s.*self.M_s).*psi_bar );
        end
        
        function psi_x = Psi_xFromNabla2PsiBar(self,nabla2_psi_bar)
            L = -sqrt(-1)*self.K./(self.K.* self.K + self.M_s.*self.M_s);
            psi_x = self.TransformBackFS( L .* nabla2_psi_bar );
        end
                
        function psi_z = Psi_zFromNabla2PsiBar(self,nabla2_psi_bar)
            L = -self.M_s./(self.K.* self.K + self.M_s.*self.M_s);
            psi_z = self.TransformBackFC( L .* nabla2_psi_bar );
        end
        
        function nabla2_psi_x = Nabla2Psi_xFromNabla2PsiBar(self,nabla2_psi_bar)
            L = sqrt(-1)*self.K;
            nabla2_psi_x = self.TransformBackFS( L .* nabla2_psi_bar );
        end
        
        function nabla2_psi_z = Nabla2Psi_zFromNabla2PsiBar(self,nabla2_psi_bar)
            L = self.M_s;
            nabla2_psi_z = self.TransformBackFC( L .* nabla2_psi_bar );
        end
        
        function b_x = b_x(self,b)
            b_x = DiffFourier(self.x,b,1,1);
        end
        
        function b_z = b_z(self,b)
            b_z = DiffSine(self.z,b,1,2);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Differential operators---Damping
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function damp_psi = damp_psi(self,nabla2_psi_bar)
            L = -(self.nu_x*self.K.^4 + self.nu_z*self.M_s.^4)./(self.K.* self.K + self.M_s.*self.M_s);
            damp_psi = self.TransformBackFS(L.*nabla2_psi_bar);
        end
        
        function [Qk,Qm] = SVV(self)
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
            
            Qk = exp( - ((abs(self.K)-k_max)./(abs(self.K)-k_cutoff)).^2 );
            Qk(abs(self.K)<k_cutoff) = 0;
            Qk(abs(self.K)>k_max) = 1;
            
            Qm = exp( - ((self.M_s-m_max)./(self.M_s-m_cutoff)).^2 );
            Qm(self.M_s<m_cutoff) = 0;
            Qm(self.M_s>m_cutoff) = 1;
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


