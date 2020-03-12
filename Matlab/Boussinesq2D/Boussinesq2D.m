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
        
        K2       
        version = 1.0
    end
        
    properties (Constant)
        g = 9.81;
    end
    
    methods(Abstract, Access = public)
        N2 = N2AtDepth(self,z)
        rho = RhoBarAtDepth(self,z)
    end
    
    methods (Abstract)%, Access = protected)
        F = InternalUVModeAtDepth(self, z, iMode) % Returns normal modes at requested depth, size(F) = [length(z) nIntModes]
        G = InternalWModeAtDepth(self, z, iMode) % Returns normal modes at requested depth, size(G) = [length(z) nIntModes]
%         F = ExternalUVModeAtDepth(self, z, iMode) % Returns normal mode at requested depth
%         G = ExternalWModeAtDepth(self, z, iMode) % Returns normal mode at requested depth 
        u = TransformToSpatialDomainWithF(self, u_bar) % Transform from (k,l,j) to (x,y,z)
        w = TransformToSpatialDomainWithG(self, w_bar ) % Transform from (k,l,j) to (x,y,z)
        ratio = UmaxGNormRatioForWave(self,k0, l0, j0) % Return the ratio/scaling required to convert a mode from the G_norm to the U_max norm
    end
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = InternalWaveModel(dims, n, z, N2)
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
            dm = 2*pi/(2*(self.Nz-1)*(t(2)-t(1)));
            self.m = dm*(0:(self.Nz-1))';
            self.m_s = self.m(2:end-1);

            [self.K,self.M] = ndgrid(self.k,self.m);
            [self.K,self.M_s] = ndgrid(self.k,self.m_s);
            [self.X,self.Z] = ndgrid(self.x,self.z);
                                    
            self.K2 = self.K.*self.K + self.M.*self.M;   % Square of the horizontal wavenumber
            self.Kh = sqrt(self.K2);
        end
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Create a single wave (public)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function period = InitializeWithPlaneWave(self, k0, j0, UAmp, sign)
            omega = self.SetGriddedWavesWithWavemodes(k0,j0,0,UAmp,sign);
            period = 2*pi/abs(omega);
        end
        
        function [f_psi, f_rho] = fluxFrom(nabla2_psi,rho)
            
        end
        
        function psi_bar = TransformForwardFS(self,psi)
            % Transform from physical coordinates to a (fourier,sine) basis
            psi_bar = SineTransformForward(self.z,psi,2,'both');
            psi_bar = FourierTransformForward(self.x,psi_bar,1);
        end
        
        function psi = TransformBackFS(self,psi_bar)
            % Transform from (fourier,sine) basis to physical coordinates
            psi_bar = SineTransformBack( self.m_s, psi_bar, 2 );
            psi = FourierTransformBack( self.k, psi_bar, 1 );
        end
        
        function psi = TransformBackFC(self,psi_bar)
            % Transform from (fourier,cosine) basis to physical coordinates
            psi_bar = cat(2,zeros(self.Nx,1),psi_bar,zeros(self.Nx,1));
            psi_bar = CosineTransformBack( self.m, psi_bar, 2);
            psi = FourierTransformBack( self.k, psi_bar, 1 );
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


