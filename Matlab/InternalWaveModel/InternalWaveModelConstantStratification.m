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
    end
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = InternalWaveModelConstantStratification(dims, n, latitude, N0)
            if length(dims) ~=3 || length(n) ~= 3
                error('The dims and n variables must be of length 3. You need to specify x,y,z');
            end
            
            if mod(log2(n(3)),1) == 0
                nz = n(3);
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
            
            N2 = N0*N0*ones(size(z));
            
            self@InternalWaveModel(dims, n, z, N2, nModes, latitude);
            
            self.N0 = N0;
            self.nz = nz;
            
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
    end
    
    methods (Access = protected)
        
        function ratio = UmaxGNormRatioForWave(self,k0, l0, j0)
            myH = self.h(k0+1,l0+1,j0);
            m = j0*pi/self.Lz;
            g = 9.81;
            F_coefficient = myH * m * sqrt(2*g/self.Lz)/sqrt(self.N0^2 - self.f0^2);
            ratio = sqrt(myH)/F_coefficient;
        end     
        
        function [F,G,h] = ModesAtWavenumber(self, k, norm)
            g = 9.81;
            k_z = (1:self.nModes)*pi/self.Lz;
            h = (self.N0*self.N0 - self.f0*self.f0)./(g*(k*k+k_z.*k_z));
            [F,G] = self.ConstantStratificationModesWithEigenvalue(k_z,h, norm);
        end
        
        function [F,G,h] = ModesAtFrequency(self, omega, norm)
            g = 9.81;
            k_z = (1:self.nModes)*pi/self.Lz;
            h = (self.N0*self.N0 - omega.*omega)./(g * k_z.*k_z);
            [F,G] = self.ConstantStratificationModesWithEigenvalue(k_z,h, norm);
        end
        
        % k_z and h should be of size [1, nModes]
        % [F,G] will return with size [length(z), nModes]
        function [F,G] = ConstantStratificationModesWithEigenvalue(self, k_z, h, norm)
            g = 9.81;
            if strcmp(norm, 'const_G_norm')
                G = sqrt(2*g/(self.Lz*(self.N0*self.N0-self.f0*self.f0))) * sin(k_z .* self.z);
                F = sqrt(2*g/(self.Lz*(self.N0*self.N0-self.f0*self.f0))) * repmat(h.*k_z,length(self.z),1) .* cos(k_z .* self.z);
            elseif strcmp(norm, 'const_F_norm')
                G = sqrt(2) * sin(k_z.*self.z) ./ repmat(h.*k_z,length(self.z),1);
                F = sqrt(2) * cos(k_z.*self.z);
            elseif strcmp(norm, 'max_u')
                G = sin(k_z.*self.z) ./ repmat(h.*k_z,length(self.z),1);
                F = cos(k_z.*self.z);
            end
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
            u = u(:,:,1:obj.Nz); % Here we use Nz (not nz) because the user may want the end point.
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
            w = w(:,:,1:obj.Nz); % Here we use Nz (not nz) because the user may want the end point.
        end
        

    end
end


