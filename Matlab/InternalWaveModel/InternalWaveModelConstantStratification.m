classdef InternalWaveModelConstantStratification < InternalWaveModel
    % InternalWaveModelConstantStratification This implements a simple
    % internal wave model for constant stratification.
    %
    % The usage is simple. First call,
    %   wavemodel = InternalWaveModelConstantStratification(dims, n, latitude, N0);
    % to initialize the model with,
    %   dims        a vector containing the length scales of x,y,z
    %   n           a vector containing the number of grid points of x,y,z
    %   latitude    the latitude of the model (e.g., 45)
    %   N0          the buoyancy frequency of the stratification
    %
    % You must now intialize the model by calling either,
    %   wavemodel.InitializeWithPlaneWave(k0, l0, j0, UAmp, sign);
    % or
    %   wavemodel.InitializeWithGMSpectrum(Amp);
    % where Amp sets the relative GM amplitude.
    %
    % Finally, you can compute u,v,w,zeta at time t by calling,
    %   [u,v] = wavemodel.VelocityFieldAtTime(t);
    %   [w,zeta] = wavemodel.VerticalFieldsAtTime(t);
    %
    % The vertical dimension must have Nz = 2^n or Nz = 2^n + 1 points. If
    % you request the extra point, then the upper boundary will be returned
    % as well. This is designed to match the DCT used by Kraig Winters'
    % model.
    %
    %   See also INTERNALWAVEMODEL and
    %   INTERNALWAVEMODELARBITRARYSTRATIFICATION
    %
    % Jeffrey J. Early
    % jeffrey@jeffreyearly.com
    properties (Access = public)
        N0
        F, G, M
    end
    
    properties (Access = protected)
        dctScratch, dstScratch;
        nz % DCT length in the vertical. This doesn't change if the user requests the value at the surface, but Nz will.
        F_cos_ext, G_sin_ext
    end
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = InternalWaveModelConstantStratification(dims, n, latitude, N0, rho0)    % MAS 1/11/18
            % rho0 is optional.
            if length(dims) ~=3 || length(n) ~= 3
                error('The dims and n variables must be of length 3. You need to specify x,y,z');
            end
            
            if mod(log2(n(3)),1) == 0
                nz = n(3);
                error('You are implicitly asking for periodic boundary conditions in the vertical. This is not supported.');
            elseif mod(log2(n(3)-1),1) == 0 % user wants the surface point
                nz = n(3)-1; % internally we proceed as if there are n-1 points
            else
                error('The vertical dimension must have 2^n or (2^n)+1 points. This is an artificial restriction.');
            end
            
            % Construct the vertical dimension
            Lz = dims(3);
            Nz = n(3);
            dz = Lz/nz;
            z = dz*(0:Nz-1)' - Lz; % cosine basis (not your usual dct basis, howeve
            
            % Number of modes is fixed for this model.
            nModes = nz;
            
            self@InternalWaveModel(dims, n, z, N0*N0*ones(size(z)), nModes, latitude);
            
            if exist('rho0','var')
                self.rho0 = rho0;
            else
                self.rho0 = 1025;
            end
            
            self.N0 = N0;
            self.nz = nz;
            
            rhoFunction = @(z) -(self.N0*self.N0*self.rho0/9.81)*z + self.rho0;
            self.internalModes = InternalModesConstantStratification([N0 self.rho0], [-dims(3) 0],z,latitude);
            
            % Preallocate this array for a faster dct
            self.dctScratch = zeros(self.Nx,self.Ny,2*self.nz);
            self.dstScratch = complex(zeros(self.Nx,self.Ny,2*self.nz));
            
            g = 9.81;
            self.M = self.J*pi/self.Lz;        % Vertical wavenumber
            h = (1/g)*(self.N0*self.N0-self.f0*self.f0)./(self.M.*self.M+self.K2);
            self.SetOmegaFromEigendepths(h);
            
            % F contains the coefficients for the U-V modes
            self.F = (self.h.*self.M)*sqrt(2*g/(self.Lz*(self.N0*self.N0-self.f0*self.f0)));
            
            % G contains the coefficients for the W-modes
            self.G = sqrt(2*g/(self.Lz*(self.N0*self.N0-self.f0*self.f0)));
            
            signNorm = -2*(mod(self.J,2) == 1)+1;
            self.F = signNorm .* self.F;
            self.G = signNorm * self.G;
        end
        
        function [u,v,w] = AnalyticalSolutionAtFrequency(self, t, x, y, z, omega, alpha, j0, phi, U)
            m = j0*pi/self.Lz;
            gh = m.*m * (self.N0*self.N0 - omega.*omega);
            K = sqrt( (omega.*omega - self.f0*self.f0)./gh );
            k0 = K.*cos(alpha);
            l0 = K.*sin(alpha);
            
            [u,v,w] = self.AnalyticalSolution(t, x, y, z, k0, l0, m, K, alpha, omega, phi, U);
        end
        
        function u = AnalyticalSolutionVectorAtWavenumber(self, t, x, k0, l0, j0, phi, U)
            xsize = size(x);
            if xsize(1) == 3
                [u,v,w] = self.AnalyticalSolutionAtWavenumber(t, x(1,:), x(2,:), x(3,:), k0, l0, j0, phi, U);
                u = cat(1,u,v,w);
            else
                [u,v,w] = self.AnalyticalSolutionAtWavenumber(t, x(:,1), x(:,2), x(:,3), k0, l0, j0, phi, U);
                u = cat(2,u,v,w);
            end
        end
        
        function [u,v,w] = AnalyticalSolutionAtWavenumber(self, t, x, y, z, k0, l0, j0, phi, U)
            alpha=atan2(l0,k0);
            K = sqrt( k0^2 + l0^2);
            m = j0*pi/self.Lz;
            gh = (self.N0 * self.N0 - self.f0*self.f0)./( m.*m + K.*K);
            omega = sqrt( gh .* K.*K + self.f0*self.f0 );
            
            [u,v,w] = self.AnalyticalSolution(t, x, y, z, k0, l0, m, K, alpha, omega, phi, U);
        end
        
        function [u,v,w] = AnalyticalSolution(self, t, x, y, z, k0, l0, m, K, alpha, omega, phi, U)
            theta = k0.*x + l0.*y + omega.*t + phi;
            cos_theta = cos(theta);
            sin_theta = sin(theta);
            u = U*(cos(alpha)*cos_theta + (self.f0/omega)*sin(alpha)*sin_theta).*cos(m*z);
            v = U*(sin(alpha)*cos_theta - (self.f0/omega)*cos(alpha)*sin_theta).*cos(m*z);
            w = (U*K/m) * sin_theta .* sin(m*z);
        end
        
        function InitializeWithFieldsAtTime(self, t, u, v, zeta)
            % This function can be used as a wave-vortex decomposition. It
            % will *exactly* recover amplitudes being used the generate the
            % dynamical fields. For the moment I assume assuming no
            % buoyancy perturbation at the boundaries.
            ubar = self.TransformFromSpatialDomainWithF( u )./self.F;
            vbar = self.TransformFromSpatialDomainWithF( v )./self.F;
            etabar = self.TransformFromSpatialDomainWithG( zeta )./self.G;
            
            delta = sqrt(self.h).*(self.K .* ubar + self.L .* vbar)./self.Kh;
            zeta = sqrt(self.h).*(self.K .* vbar - self.L .* ubar)./self.Kh;
            
            A_plus = exp(-sqrt(-1)*self.Omega*t).*(-self.g*self.Kh.*sqrt(self.h).*etabar./self.Omega + delta - sqrt(-1)*zeta*self.f0./self.Omega)/2;
            A_minus = exp(sqrt(-1)*self.Omega*t).*(self.g*self.Kh.*sqrt(self.h).*etabar./self.Omega + delta + sqrt(-1)*zeta*self.f0./self.Omega)/2;
            B = (etabar*self.f0 - sqrt(-1)*zeta.*self.Kh.*sqrt(self.h))*self.f0./(self.Omega.*self.Omega);
            
            % inertial must be solved for separately.
            A_plus(1,1,:) = exp(-sqrt(-1)*self.f0*t)*(ubar(1,1,:) - sqrt(-1)*vbar(1,1,:)).*sqrt(self.h(1,1,:))/2;
            A_minus(1,1,:) = conj(A_plus(1,1,:));
            B(1,1,:) = 0;
            
            % B is the geostrophic solution, not yet implemented.
            A_plus = InternalWaveModel.MakeHermitian(A_plus);
            A_minus = InternalWaveModel.MakeHermitian(A_minus);
            self.GenerateWavePhases(A_plus,A_minus);
        end
        
        function InitializeWithIsopycnalDisplacementField(self, zeta)
            % Note that you will lose any 'mean' zeta, because technically
            % this is the perturbation *from* a mean.
            zeta_bar = self.TransformFromSpatialDomainWithG(zeta);
            
            Kh = self.Kh;
            Kh(Kh<1e-14) = 1;
            negate_zeta = abs(self.Omega)./(Kh .* sqrt(self.h));
            A_plus = -zeta_bar .* negate_zeta ./ self.G; % extra factor from constant stratification case.
            self.GenerateWavePhases(A_plus,zeros(size(self.K)));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computes the phase information given the amplitudes (internal)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function GenerateWavePhases(self, U_plus, U_minus)
            GenerateWavePhases@InternalWaveModel(self, U_plus, U_minus);
            
            % Multiply by the norm coefficients since we're using a
            % discrete cosine transform directly.
            self.u_plus = self.u_plus .* self.F;
            self.u_minus = self.u_minus .* self.F;
            
            self.v_plus = self.v_plus .* self.F;
            self.v_minus = self.v_minus .* self.F;
            
            self.w_plus = self.w_plus .* self.G;
            self.w_minus = self.w_minus .* self.G;
            
            self.zeta_plus = self.zeta_plus .* self.G;
            self.zeta_minus = self.zeta_minus .* self.G;
        end
                
        function PrecomputeExternalWaveCoefficients(self)
            PrecomputeExternalWaveCoefficients@InternalWaveModel(self)
            if self.norm_ext == Normalization.kConstant
                g = 9.81;
                coeff = sqrt(2*g/(self.Lz*(self.N0*self.N0-self.f0*self.f0)));
                self.F_cos_ext = coeff * ( self.h_ext .* self.k_z_ext );
                self.G_sin_ext = coeff * ones(size(self.h_ext));
            elseif self.norm_ext == Normalization.uMax
                self.F_cos_ext = ones(size(self.h_ext));
                self.G_sin_ext = 1./(self.h_ext .* self.k_z_ext);
            end
        end
        
        function rho = RhoBarAtDepth(self,z)
            g = 9.81;
            rho = -(self.N0*self.N0*self.rho0/g)*z + self.rho0;
        end
        
        function N2 = N2AtDepth(self,z)
            N2 = self.N0 * self.N0 * ones(size(z));
        end
    end
    
    methods %(Access = protected)
        
        function ratio = UmaxGNormRatioForWave(self,k0, l0, j0)
            myH = self.h(k0+1,l0+1,j0);
            m = j0*pi/self.Lz;
            g = 9.81;
            F_coefficient = myH * m * sqrt(2*g/self.Lz)/sqrt(self.N0^2 - self.f0^2);
            ratio = sqrt(myH)/F_coefficient;
        end     
                
        % size(z) = [N 1]
        % size(waveIndices) = [1 M]
        function F = ExternalUVModeAtDepth(self, z, iMode)
            % Called by the superclass when advecting particles spectrally.
            F = self.F_cos_ext(iMode) * cos( z * self.k_z_ext(iMode) );
        end
        
        function G = ExternalWModeAtDepth(self, z, iMode)
            % Called by the superclass when advecting particles spectrally.
            G = self.G_sin_ext(iMode) * sin( z * self.k_z_ext(iMode) );
        end
                
        function F = InternalUVModeAtDepth(self, z, iMode)
            % Called by the superclass when advecting particles spectrally.
            k_z = self.j_int(iMode)*pi/self.Lz;
            g = 9.81;
            coeff = sqrt(2*g/(self.Lz*(self.N0*self.N0-self.f0*self.f0)));
            F = coeff * ( self.h_int(iMode) * k_z ) .* cos(z * k_z); % [N M]
        end
        
        function G = InternalWModeAtDepth(self, z, iMode)
            % Called by the superclass when advecting particles spectrally.
            k_z = self.j_int(iMode)*pi/self.Lz;
            g = 9.81;
            coeff = sqrt(2*g/(self.Lz*(self.N0*self.N0-self.f0*self.f0)));
            G = coeff * sin(z * k_z); % [N M]
        end
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computes the phase information given the amplitudes (internal)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function u = TransformToSpatialDomainWithF(obj, u_bar)
            % Here we use what I call the 'Fourier series' definition of the ifft, so
            % that the coefficients in frequency space have the same units in time.
            u = obj.Nx*obj.Ny*ifft(ifft(u_bar,obj.Nx,1),obj.Ny,2,'symmetric');
            
            % Re-order to convert to an fast cosine transform
            obj.dctScratch = cat(3, zeros(obj.Nx,obj.Ny), 0.5*u(:,:,1:obj.nz-1), u(:,:,obj.nz), 0.5*u(:,:,obj.nz-1:-1:1));
  
            u = fft(obj.dctScratch,2*obj.nz,3);
            if obj.performSanityChecks == 1
                ratio = max(max(max(abs(imag(u)))))/max(max(max(abs(real(u)))));
                if ratio > 1e-6
                    fprintf('WARNING: The inverse cosine transform reports an unreasonably large imaginary part, %.2g.\n',ratio);
                end
            end
            % should not have to call real, but for some reason, with enough
            % points, it starts generating some small imaginary component.
            u = real(u(:,:,1:obj.Nz)); % Here we use Nz (not nz) because the user may want the end point.
        end
        
        function w = TransformToSpatialDomainWithG(obj, w_bar )
            % Here we use what I call the 'Fourier series' definition of the ifft, so
            % that the coefficients in frequency space have the same units in time.
            w = obj.Nx*obj.Ny*ifft(ifft(w_bar,obj.Nx,1),obj.Ny,2,'symmetric');
            
            % Re-order to convert to an fast cosine transform
            obj.dstScratch = sqrt(-1)*cat(3, zeros(obj.Nx,obj.Ny), 0.5*w(:,:,1:obj.nz-1), w(:,:,obj.nz), -0.5*w(:,:,obj.nz-1:-1:1));
            
            w = fft( obj.dstScratch,2*obj.nz,3);
            if obj.performSanityChecks == 1
                ratio = max(max(max(abs(imag(w)))))/max(max(max(abs(real(w)))));
                if ratio > 1e-6
                    fprintf('WARNING: The inverse sine transform reports an unreasonably large imaginary part, %.2g.\n',ratio);
                end
            end
            % should not have to call real, but for some reason, with enough
            % points, it starts generating some small imaginary component.
            w = real(w(:,:,1:obj.Nz)); % Here we use Nz (not nz) because the user may want the end point.
        end
                
        function u_bar = TransformFromSpatialDomainWithF(self, u)
            self.dctScratch = ifft(cat(3,u,u(:,:,self.nz:-1:2)),2*self.nz,3);
            u_bar = 2*real(self.dctScratch(:,:,2:self.nz+1)); % we *ignore* the barotropic mode, starting at 2, instead of 1 in the transform
            u_bar = fft(fft(u_bar,self.Nx,1),self.Ny,2)/self.Nx/self.Ny;
        end
        
        function w_bar = TransformFromSpatialDomainWithG(self, w)
            self.dstScratch = ifft(cat(3,w,-w(:,:,self.nz:-1:2)),2*self.nz,3);
            w_bar = 2*imag(self.dstScratch(:,:,2:self.nz+1));
            w_bar = fft(fft(w_bar,self.Nx,1),self.Ny,2)/self.Nx/self.Ny;
        end

    end
end


