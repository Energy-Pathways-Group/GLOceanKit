%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% InternalWaveModel
%
% This implements a simple internal wave model for constant stratification.
%
% The usage is simple. First call,
%   wavemodel = InternalWaveModel(dims, n, latitude, N0);
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
% The vertical dimension must have Nz = 2^n or Nz = 2^n + 1 points. If you
% request the extra point, then the upper boundary will be returned as
% well. This is designed to match the DCT used by Kraig Winters' model.
%
%
% Jeffrey J. Early
% jeffrey@jeffreyearly.com
%
% March 25th, 2016      Version 1.0
% March 30th, 2016      Version 1.1
% November 17th, 2016   Version 1.2
% December 9th, 2016    Version 1.3
% February 9th, 2017    Version 1.4

classdef InternalWaveModelConstantStratification < InternalWaveModel
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

        function ratio = UmaxGNormRatioForWave(self,k0, l0, j0)
            myH = self.h(k0+1,l0+1,j0);
            m = j0*pi/self.Lz;
            g = 9.81;
            F_coefficient = myH * m * sqrt(2*g/self.Lz)/sqrt(self.N0^2 - self.f0^2);
            ratio = sqrt(myH)/F_coefficient;
        end
        
        function self = SetExternalWavesWithWavenumbers(self, k, l, j, phi, U)
            self.kExternal = k;
            self.lExternal = l;
            self.alphaExternal = atan2(l,k);
            self.phiExternal = phi;
            self.uExternal = U;
            
            g = 9.81;
            for iWave=1:length(j)
                kz = j(iWave)*pi/self.Lz;        % Vertical wavenumber
                K2h = k(iWave)*k(iWave) + l(iWave)*l(iWave);
                c2 = (self.N0*self.N0-self.f0*self.f0)./(kz*kz + K2h);
                self.omegaExternal(iWave) = sqrt( c2 * k2h + self.f0*self.f0 );
                
            end
        end
        
        function self = SetExternalWavesWithFrequencies(self, omega, alpha, j, phi, U)
            self.omegaExternal = omega;
            self.alphaExternal = alpha;
            self.phiExternal = phi;
            self.uExternal = U;
            
            for iWave=1:length(j)
                kz = j(iWave)*pi/self.Lz;        % Vertical wavenumber  
                h = (1/g)*(self.N0*self.N0-self.f0*self.f0)./(kz*kz+self.K2);
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


