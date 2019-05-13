classdef (Abstract) InternalWaveModel < handle
    % InternalWaveModel This abstract class models linear internal waves.
    % There are two subclasses: one for arbitrary stratification profiles
    % and one for constant stratification profiles (which can therefore use
    % fast transforms in the vertical).
    %
    % The usage is simple. First call,
    %   wavemodel = InternalWaveModelConstantStratification(dims, n, latitude, N0);
    % to initialize the model with constant stratification or,
    %   wavemodel = InternalWaveModelArbitraryStratification(dims, n, rho, z, nModes, latitude)
    % to initialize the model with artibrary stratification given by rho.
    %
    % By default the model will include all possible waves allowed by an FFT
    % given the horizontal grid size. However, in some cases you may want to
    % include *additional* waves (of say, longer wavelength) than otherwise
    % would be allowed. In that case you may call
    %   wavemodel.SetExternalWavesWithWavenumbers(self, k, l, j, phi, A, type)
    % or,
    %   wavemodel.SetExternalWavesWithFrequencies(self, omega, alpha, j, phi, A, type)
    % to added these additional waves. These additional waves cannot take
    % advantage of fast transforms, and therefore adding too many will
    % significantly slow the model.
    %
    % Note that if you make the aspect ratio something other than 1, you
    % want more points in the first (x) dimension. That makes for faster
    % FFTs.
    %
    % You must now intialize the model by calling either,
    %   wavemodel.InitializeWithPlaneWave(k0, l0, j0, UAmp, sign);
    % or
    %   wavemodel.InitializeWithGMSpectrum(Amp);
    % where Amp sets the relative GM amplitude.
    %
    % Finally, you can compute u,v,w,rho at time t by calling,
    %   [u,v,w] = wavemodel.VelocityFieldAtTime(t);
    %   rho = wavemodel.DensityFieldAtTime(t);
    %
    %   See also INTERNALWAVEMODELCONSTANTSTRATIFICATION and
    %   INTERNALWAVEMODELARBITRARYSTRATIFICATION
    %
    % Jeffrey J. Early
    % jeffrey@jeffreyearly.com
    %
    % March 25th, 2016      Version 1.0
    % March 30th, 2016      Version 1.1
    % November 17th, 2016   Version 1.2
    % December 9th, 2016    Version 1.3
    % February 9th, 2017    Version 1.4
    % March 20th, 2017      Version 1.5
    % December 15th, 2017   Version 1.6
    properties (Access = public)
        Lx, Ly, Lz % Domain size
        Nx, Ny, Nz % Number of points in each direction
        nModes
        latitude
        
        x, y, z
        k, l, j
        X,Y,Z
        K,L,J
        
        rhobar, rho0
        N2, Nmax
        
        K2, Kh, h, Omega, f0, C
        u_plus, u_minus, v_plus, v_minus, w_plus, w_minus, zeta_plus, zeta_minus
        
        Amp_plus, Amp_minus % Used for diagnostics only
        B % amplitude of the geostrophic component (matches Amp_plus/minus)
        B0 % amplitude of the geostrophic component for k_z=0.

        
        didPreallocateAdvectionCoefficients = 0
        U_cos_int, U_sin_int, V_cos_int, V_sin_int, W_sin_int, Zeta_cos_int, Rho_cos_int, k_int, l_int, j_int, h_int, omega_int, phi_int
        u_g, v_g, zeta_g % geostropic components
        
        % These are all row vectors, e.g. size(U_ext)=[1 length(U_ext)], except F_ext, G_ext which are size(F_ext) = [length(z) length(U_ext)];
        U_ext, k_ext, l_ext, j_ext, k_z_ext, h_ext, omega_ext, phi_ext, F_ext, G_ext, norm_ext
        U_cos_ext, U_sin_ext, V_cos_ext, V_sin_ext, W_sin_ext, Zeta_cos_ext
        
        Xc, Yc, Zc % These may contain 'circular' versions of the grid
                
        internalModes % InternalModes object being used by the model
        
        version = 1.6
        performSanityChecks = 0
        advectionSanityCheck = 0
    end
    
    properties (Dependent)
        Rho_cos_ext % Don't remove. Cim uses this.
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
        function self = InternalWaveModel(dims, n, z, N2, nModes, latitude)
            if length(dims) ~=3 || length(n) ~= 3
                error('The dims and n variables must be of length 3. You need to specify x,y,z');
            end
            
            self.Lx = dims(1);
            self.Ly = dims(2);
            self.Lz = dims(3);
            
            self.Nx = n(1);
            self.Ny = n(2);
            self.Nz = n(3);
            self.nModes = nModes;
            
            self.latitude = latitude;
            
            dx = self.Lx/self.Nx;
            dy = self.Ly/self.Ny;
            
            self.x = dx*(0:self.Nx-1)'; % periodic basis
            self.y = dy*(0:self.Ny-1)'; % periodic basis
            self.z = z; % cosine basis (not your usual dct basis, however)
            
            self.N2 = N2;
            self.Nmax = sqrt(max(N2));
%             self.rhobar = self.RhoBarAtDepth(self.z);
            
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
            
            self.j = (1:nModes)';
            
            [K,L,J] = ndgrid(self.k,self.l,self.j);
            [X,Y,Z] = ndgrid(self.x,self.y,self.z);
                        
            self.L = L; self.K = K; self.J = J;
            self.X = X; self.Y = Y; self.Z = Z;
            
            self.f0 = 2 * 7.2921E-5 * sin( self.latitude*pi/180 );
            self.K2 = self.K.*self.K + self.L.*self.L;   % Square of the horizontal wavenumber
            self.Kh = sqrt(self.K2);
            
            nothing = zeros(size(self.K));
            self.h = nothing; self.Omega = nothing;
            self.C = nothing; self.u_plus = nothing;
            self.u_minus = nothing; self.v_plus = nothing; self.v_minus = nothing;
            self.w_plus = nothing; self.w_minus = nothing; self.zeta_plus = nothing;
            self.zeta_minus = nothing;
            self.Amp_minus = nothing; self.Amp_plus = nothing;
        end
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Create a single wave (public)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function period = InitializeWithPlaneWave(self, k0, l0, j0, UAmp, sign)
            omega = self.SetGriddedWavesWithWavemodes(k0,l0,j0,0,UAmp,sign);
            period = 2*pi/abs(omega);
        end
        
        function RemoveAllGriddedWaves(self)
            self.GenerateWavePhases(zeros(size(self.K)),zeros(size(self.K)));            
            self.didPreallocateAdvectionCoefficients = 0;
        end
        
        function [omega,k,l] = SetGriddedWavesWithWavemodes(self, kMode, lMode, jMode, phi, Amp, signs)
            self.RemoveAllGriddedWaves();
            [omega,k,l] = self.AddGriddedWavesWithWavemodes(kMode, lMode, jMode, phi, Amp, signs);
        end
        
        function [omega,k,l] = AddGriddedWavesWithWavemodes(self, kMode, lMode, jMode, phi, Amp, signs)
            % Add wavemodes on the gridded field.
            % The values given must meet the following requirements:
            % (k0 > -Nx/2 && k0 < Nx/2)
            % (l0 > -Ny/2 && l0 < Ny/2)
            % (j0 >= 1 && j0 < nModes)
            % phi is in radians, from 0-2pi
            % Amp is the fluid velocity U
            % sign is +/-1, indicating the sign of the frequency.
            
            if ~isequal(size(kMode), size(lMode), size(jMode), size(phi), size(Amp), size(signs))
                error('All input array must be of equal size');
            end
            A_minus_total = zeros(size(self.K));
            A_plus_total = zeros(size(self.K));
            omega = zeros(size(kMode));
            k = zeros(size(kMode));
            l = zeros(size(kMode));
            for iMode = 1:length(kMode)
                k0 = kMode(iMode);
                l0 = lMode(iMode);
                j0 = jMode(iMode);
                phi0 = phi(iMode);
                A = Amp(iMode);
                sign = signs(iMode);
                
                % User input sanity checks. We don't deal with the Nyquist.
                if (k0 <= -self.Nx/2 || k0 >= self.Nx/2)
                    error('Invalid choice for k0 (%d). Must be an integer %d < k0 < %d',k0,-self.Nx/2+1,self.Nx/2-1);
                end
                if (l0 <= -self.Ny/2 || l0 >= self.Ny/2)
                    error('Invalid choice for l0 (%d). Must be an integer %d < l0 < %d',l0,-self.Ny/2+1,self.Ny/2+1);
                end
                if (j0 < 1 || j0 >= self.nModes)
                    error('Invalid choice for j0 (%d). Must be an integer 0 < j < %d',j0, self.nModes);
                end
                
                % Deal with the negative wavenumber cases (and inertial)
                if l0 == 0 && k0 == 0 % inertial
                    if sign < 1
                        sign=1;   
                        phi0 = -phi0;
                    end
                elseif l0 == 0 && k0 < 0
                    k0 = -k0;
                    sign = -1*sign;
                    A = -1*A;
                    phi0 = -phi0;
                elseif l0 < 0
                    l0 = -l0;
                    k0 = -k0;
                    sign = -1*sign;
                    A = -1*A;
                    phi0 = -phi0;
                end
                
                % Rewrap (k0,l0) to follow standard FFT wrapping. l0 should
                % already be correct.
                if (k0 < 0)
                    k0 = self.Nx + k0;
                end
                
                ratio = self.UmaxGNormRatioForWave(k0, l0, j0);
                
                U = zeros(size(self.K));
                U(k0+1,l0+1,j0) = A*ratio/2*exp(sqrt(-1)*phi0);
                if sign > 0
                    A_plus = InternalWaveModel.MakeHermitian(U);
                    A_minus = zeros(size(U));
                    A_minus(1,1,:) = conj(A_plus(1,1,:)); % Inertial oscillations are created using this trick.
                else
                    A_plus = zeros(size(U));
                    A_minus = InternalWaveModel.MakeHermitian(U);
                end
                
                A_minus_total = A_minus_total + A_minus;
                A_plus_total = A_plus_total + A_plus;
                
                % When we hand back the actual frequency and wavenumbers,
                % we honor the users original intent and the match the
                % signs they provided.
                if (kMode(iMode) < 0)
                    k_out = self.Nx + kMode(iMode);
                else
                    k_out = kMode(iMode);
                end
                if (lMode(iMode) < 0)
                    l_out = self.Ny + lMode(iMode);
                else
                    l_out = lMode(iMode);
                end
                omega(iMode) = signs(iMode)*abs(self.Omega(k0+1,l0+1,j0));
                k(iMode) = self.K(k_out+1,l_out+1,j0);
                l(iMode) = self.L(k_out+1,l_out+1,j0);
            end
            
            % We literally just add this wave to the existing waves.
            self.GenerateWavePhases(self.Amp_plus + A_plus_total,self.Amp_minus + A_minus_total);
            
            self.didPreallocateAdvectionCoefficients = 0;
        end
        
        function [omega, alpha, k, l, mode, phi, A, norm] = WaveCoefficientsFromGriddedWaves(self)

            % This returns the properties of the waves being used in the
            % gridded simulation, as their properly normalized individual
            % wave components. Very useful for debugging.
            %
            % Note that A_plus and A_minus each have half the inertial
            % energy. This can be misleading, but the phasing is chosen to
            % make it work. Never-the-less, we double/zero that component.
            A_p = self.Amp_plus;
            A_p(1,1,:) = 2*A_p(1,1,:);
            A_m = self.Amp_minus;
            A_m(1,1,:) = 0*A_m(1,1,:);
            
            [A_plus,phi_plus,linearIndex] = ExtractNonzeroWaveProperties(A_p);
            omega_plus = self.Omega(linearIndex);
            mode_plus = self.J(linearIndex);
            alpha_plus = atan2(self.L(linearIndex),self.K(linearIndex));
            k_plus = self.K(linearIndex);
            l_plus = self.L(linearIndex);
            
            [A_minus,phi_minus,linearIndex] = ExtractNonzeroWaveProperties(A_m);
            omega_minus = -self.Omega(linearIndex);
            mode_minus = self.J(linearIndex);
            alpha_minus = atan2(self.L(linearIndex),self.K(linearIndex));
            k_minus = self.K(linearIndex);
            l_minus = self.L(linearIndex);
            
            k = [k_plus; k_minus];
            l = [l_plus; l_minus];
            omega = [omega_plus; omega_minus];
            mode = [mode_plus; mode_minus];
            alpha = [alpha_plus; alpha_minus];
            phi = [phi_plus; phi_minus];
            A = [A_plus; A_minus];
            norm = Normalization.kConstant;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % These functions add free waves to the model, not constrained to
        % the FFT grid.
        %
        % Add free waves to the model that are not constrained to the
        % gridded solution. The amplitude, A, can be given as an energy
        % density or a maximum U velocity. The type must then be either
        % Normalization.kConstant (energy density) or Normalization.uMax.
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function RemoveAllExternalWaves(self)
            % Remove all external waves.
            self.U_ext = [];
            self.k_ext = [];
            self.l_ext = [];
            self.j_ext = [];
            self.k_z_ext = [];
            self.phi_ext = [];
            self.h_ext = [];
            self.F_ext = [];
            self.G_ext = [];
            self.omega_ext = [];
            self.PrecomputeExternalWaveCoefficients(); 
        end
        
        function omega = SetExternalWavesWithWavenumbers(self, k, l, j, phi, A, norm)
            % Replaces the existing set of external modes with new ones
            % given with horizontal wavenumbers (k,l), in radians per
            % meter.
            %
            % j indicates the vertical mode number (j>=1) phi indicates the
            % phase of the wave, in radians A indicates the amplitude of
            % the wave, with respect to the given norm, which should be
            % either Normalization.uMax or Normalization.kConstant.
            self.RemoveAllExternalWaves();
            omega = self.AddExternalWavesWithWavenumbers(k, l, j, phi, A, norm);
        end
        
        function omega = AddExternalWavesWithWavenumbers(self, k, l, j, phi, A, norm)
            % Adds external modes with horizontal wavenumbers (k,l), in
            % radians per meter.
            %
            % j indicates the vertical mode number (j>=1)
            % phi indicates the phase of the wave, in radians
            % A indicates the amplitude of the wave, with respect to the
            % given norm, which should be either Normalization.uMax or
            % Normalization.kConstant.
            if ~isequal(size(k), size(l), size(j), size(phi), size(A))
                error('All input array must be of equal size');
            end
            K2h = reshape(k.*k + l.*l,1,[]);  
            [h_, validIndices] = self.AddExternalWavesWithMethod(j,phi,A,norm,sqrt(K2h),'ModesAtWavenumber');
            K2h = K2h(validIndices);
            k = k(validIndices);
            l = l(validIndices);
            
            if ~isempty(k)
                omega = sqrt(self.g*h_ .* K2h + self.f0*self.f0);
                
                self.k_ext = cat(2,self.k_ext,reshape(k,1,[]));
                self.l_ext = cat(2,self.l_ext,reshape(l,1,[]));
                self.omega_ext = cat(2,self.omega_ext,reshape(omega,1,[]));
                
                self.PrecomputeExternalWaveCoefficients();
            end
        end
        
        function k = SetExternalWavesWithFrequencies(self, omega, alpha, j, phi, A, norm)
            % Replaces the existing set of external modes with new ones
            % given with frequency omega (radians/second) and phase angle
            % alpha (radians).
            %
            % j indicates the vertical mode number (j>=1)
            % phi indicates the phase of the wave, in radians
            % A indicates the amplitude of the wave, with respect to the
            % given norm, which should be either Normalization.uMax or
            % Normalization.kConstant.
            self.RemoveAllExternalWaves();
            k = self.AddExternalWavesWithFrequencies(omega, alpha, j, phi, A, norm);
        end
        
        function k = AddExternalWavesWithFrequencies(self, omega, alpha, j, phi, A, norm)
            % Adds external modes with frequency omega (radians/second) and
            % phase angle alpha (radians).
            %
            % j indicates the vertical mode number (j>=1)
            % phi indicates the phase of the wave, in radians
            % A indicates the amplitude of the wave, with respect to the
            % given norm, which should be either Normalization.uMax or
            % Normalization.kConstant.
            if ~isequal(size(omega), size(alpha), size(j), size(phi), size(A))
                error('All input array must be of equal size');
            end
            omega = reshape(omega,1,[]);
            [h_, validIndices] = self.AddExternalWavesWithMethod(j,phi,A,norm,omega,'ModesAtFrequency');
            
            omega = omega(validIndices);
            alpha = alpha(validIndices);
            
            if ~isempty(omega)
                k = sqrt((omega.*omega - self.f0*self.f0)./(self.g*h_));
                alpha0 = reshape(alpha,1,[]);
                
                self.k_ext = cat(2,self.k_ext,reshape(k .* cos(alpha0),1,[]));
                self.l_ext = cat(2,self.l_ext,reshape(k .* sin(alpha0),1,[]));
                self.omega_ext = cat(2,self.omega_ext,omega);
                
                self.PrecomputeExternalWaveCoefficients();
            else
                k = [];
            end
        end
        
        
        function [h, validIndices] = AddExternalWavesWithMethod( self, j, phi, A, norm, kOrOmega, methodName )
            % This function is called by AddExternalWavesWithFrequencies
            % and AddExternalWavesWithWavenumbers and should not be used
            % directly.
            %
            % The function returns a culled list of h thats contain only
            % valid values, and a list of the validIndices from the
            % original set.
            %
            % Appends to j_ext, k_z_ext, phi_ext, F_ext, G_ext, h_ext, and
            % U_ext. The calling functions are still responsible for
            % setting k_ext, l_ext, and omega_ext
            numExistingWaves = length(self.k_ext);
            validIndices = zeros(size(kOrOmega));
            h = zeros(size(kOrOmega));
            
            switch norm
                case {Normalization.kConstant, Normalization.uMax}
                    self.norm_ext = norm;
                otherwise
                    error('Invalid norm. You must use Normalization.kConstant or Normalization.uMax.');
            end
            
            self.internalModes.normalization = self.norm_ext;
            numValidIndices = 0;
            for iWave=1:length(kOrOmega)
                [FExt,GExt,hExt] = self.internalModes.(methodName)(abs(kOrOmega(iWave)));
                if (hExt(j(iWave)) <= 0)
                    warning('You attempted to add a wave that returned an invalid eigenvalue! It will be skipped. You tried to add the j=%d mode computed with %s=%f which returned eigenvalue h=%f.\n', j(iWave), methodName, kOrOmega(iWave), hExt(j(iWave)));
                    continue;
                end
                numValidIndices = numValidIndices + 1;
                validIndices(numValidIndices) = iWave;
                h(numValidIndices) = hExt(j(iWave));
                
                index = numExistingWaves + numValidIndices;
                                
                self.F_ext(:,index) = FExt(:,j(iWave));
                self.G_ext(:,index) = GExt(:,j(iWave));
                self.h_ext = cat(2,self.h_ext,hExt(j(iWave)));
                self.j_ext = cat(2,self.j_ext,j(iWave));
                self.k_z_ext = cat(2,self.k_z_ext,j(iWave)*pi/self.Lz);
                self.phi_ext = cat(2,self.phi_ext,phi(iWave));
                if norm == Normalization.kConstant
                    self.U_ext = cat(2,self.U_ext,A(iWave)/sqrt(hExt(j(iWave))));
                elseif norm == Normalization.uMax
                    self.U_ext = cat(2,self.U_ext,A(iWave));
                end
            end
            
            validIndices = validIndices(1:numValidIndices);
            h = h(1:numValidIndices);
        end
        
        function PrecomputeExternalWaveCoefficients(self)
            alpha0 = atan2(self.l_ext,self.k_ext);
            Kh_ = sqrt( self.k_ext.*self.k_ext + self.l_ext.*self.l_ext);
            
            kOverOmega = Kh_ ./ self.omega_ext;
            if self.latitude == 0
                f0OverOmega = 0;
                kOverOmega( Kh_ == 0 ) = 0;
            else
                f0OverOmega = (self.f0 ./ self.omega_ext);  
            end
            
            self.U_cos_ext = self.U_ext .* cos(alpha0);
            self.U_sin_ext = self.U_ext .* f0OverOmega .* sin(alpha0);
            self.V_cos_ext = self.U_ext .* sin(alpha0);
            self.V_sin_ext = -self.U_ext .* f0OverOmega .* cos(alpha0);
            self.W_sin_ext = self.U_ext .* Kh_ .* self.h_ext;
            self.Zeta_cos_ext = - self.U_ext .* kOverOmega .* self.h_ext;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Create a full Garrett-Munk spectrum (public)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function InitializeWithGMSpectrum(self, GMAmplitude, varargin)
            if mod(length(varargin),2) ~= 0
                error('Arguments must be given as name/value pairs.');
            end
            
            % Set defaults
            j_star = 3;
            
            % Now override the defaults with user settings
            for iArg = 1:2:length(varargin)
                if strcmp(varargin{iArg}, 'j_star')
                    j_star = varargin{iArg+1};
                    varargin(iArg+1) = [];
                    varargin(iArg) = [];
                    break;
                end
            end
            
            % GM Parameters
            L_gm = 1.3e3; % thermocline exponential scale, meters
            invT_gm = 5.2e-3; % reference buoyancy frequency, radians/seconds
            E_gm = 6.3e-5; % non-dimensional energy parameter
            E = L_gm*L_gm*L_gm*invT_gm*invT_gm*E_gm*GMAmplitude;
%             E = E*(self.Lz/L_gm); % This correction fixes the amplitude so that the HKE variance at a given depth matches (instead of depth integrated energy)
                      
            % Compute the proper vertical function normalization
            H = (j_star+(1:1024)).^(-5/2);
            H_norm = 1/sum(H);
            
            % Do the same for the frequency function.
            B_norm = 1/atan(sqrt(self.Nmax*self.Nmax/(self.f0*self.f0)-1));
            
            % This function tells you how much energy you need between two
            % frequencies for a given vertical mode.
            GM2D_int = @(omega0,omega1,j) E*H_norm*B_norm*((j+j_star).^(-5/2))*(atan(self.f0/sqrt(omega0*omega0-self.f0*self.f0)) - atan(self.f0/sqrt(omega1*omega1-self.f0*self.f0)));
            
            % Do a quick check to see how much energy is lost due to
            % limited vertical resolution.
            maxMode = self.nModes;
            for iArg = 1:2:length(varargin)
                if strcmp(varargin{iArg}, 'maxMode')
                    maxMode = varargin{iArg+1};
                end
            end
                  
            totalEnergy = 0;
            for mode=1:maxMode
                totalEnergy = totalEnergy + GM2D_int(self.f0,self.Nmax,mode);
            end
            fprintf('You will miss %.2f%% of the energy due to limited vertical modes.\n',100-100*totalEnergy/E);
            
            [GM3Dint,GM3Dext] = self.InitializeWithSpectralFunction(GM2D_int,varargin{:});
            
            fprintf('After distributing energy across frequency and mode, you still have %.2f%% of reference GM energy.\n',100*(sum(sum(sum(GM3Dint))) + sum(GM3Dext))/E);
            fprintf('Due to restricted domain size, the j=1,k=l=0 mode contains %.2f%% the total energy.\n',100*GM3Dint(1,1,1)/(sum(sum(sum(GM3Dint))) + sum(GM3Dext)) );
            
            GM_sum_int = sum(sum(sum(GM3Dint)))/E;
            GM_sum_ext = sum(GM3Dext)/E;
            GM_random_sum_int = sum(sum(sum(self.Amp_plus.*conj(self.Amp_plus) + self.Amp_minus.*conj(self.Amp_minus)  )))/E;
            GM_random_sum_ext = sum((self.U_cos_ext.*self.U_cos_ext + self.V_cos_ext.*self.V_cos_ext).*self.h_ext/2)/E;
            fprintf('The (gridded, external) wave field sums to (%.2f%%, %.2f%%) GM given the scales, and the randomized field sums to (%.2f%%, %.2f%%) GM\n', 100*GM_sum_int, 100*GM_sum_ext, 100*GM_random_sum_int,100*GM_random_sum_ext);
        end
        
        function [GM3Dint,GM3Dext] = InitializeWithSpectralFunction(self, GM2D_int, varargin)   
            % The GM2D_int function is used to assign variance to a given
            % wave mode. It has three arguments, omega0, omega1, and j and
            % should return the amount of variance you want assigned to a
            % wave mode between omega0 and omega1 at vertical mode j.
            %
            % The returned values GM3Dint are the results of distributing
            % this variance. size(GM3Dint) = size(self.Kh), so you can see
            % how much energy was assigned to each internal mode and
            % similarly size(GM3Dext) = size(self.k_ext).
            %
            % The function takes the (optional) name/value pairs:
            %
            % shouldRandomizeAmplitude = 1 or 0 will randomize the
            % energy in each mode such that the expected value matches that
            % assigned. Default 0 (amplitudes will not be randomized)
            %
            % maxDeltaOmega is the maximum width in frequency that will be
            % integrated over for assigned energy. By default it is self.Nmax-self.f0
            if nargin(GM2D_int) ~= 3
                error('The spectral function must take three inputs: omega0, omega1, and j.\n');
            end
            
            if mod(length(varargin),2) ~= 0
                error('Arguments must be given as name/value pairs.');
            end
            
            % Set defaults
            shouldRandomizeAmplitude = 0;
            maxDeltaOmega = self.Nmax-self.f0;
            initializeModes = 0;
            energyWarningThreshold = 0.5;
            excludeNyquist = 1;
            minK = 0;
            maxK = max(max(max(abs(self.K))));
            minMode = 1;
            maxMode = self.nModes;
            
            % Now override the defaults with user settings
            for iArg = 1:2:length(varargin)
                if strcmp(varargin{iArg}, 'shouldRandomizeAmplitude')
                    shouldRandomizeAmplitude = varargin{iArg+1};
                elseif strcmp(varargin{iArg}, 'maxDeltaOmega')
                    maxDeltaOmega = varargin{iArg+1};
                elseif strcmp(varargin{iArg}, 'energyWarningThreshold')
                    energyWarningThreshold = varargin{iArg+1};
                elseif strcmp(varargin{iArg}, 'excludeNyquist')
                    excludeNyquist = varargin{iArg+1};
                elseif strcmp(varargin{iArg}, 'initializeModes')
                    if strcmp(varargin{iArg+1}, 'all')
                        initializeModes = 0;
                    elseif strcmp(varargin{iArg+1}, 'internalOnly')
                        initializeModes = 1;
                    elseif strcmp(varargin{iArg+1}, 'externalOnly')
                        initializeModes = 2;
                    else
                        error('Invalid option for initializeModes');
                    end
                elseif strcmp(varargin{iArg}, 'minK')
                    minK = varargin{iArg+1};
                elseif strcmp(varargin{iArg}, 'maxK')
                    maxK = varargin{iArg+1};
                elseif strcmp(varargin{iArg}, 'minMode')
                    minMode = varargin{iArg+1};
                elseif strcmp(varargin{iArg}, 'maxMode')
                    maxMode = varargin{iArg+1};
                else
                    error('Invalid argument');
                end
            end
                  
            if excludeNyquist == 1
                nyquistIndicesForK = sub2ind(size(self.Omega),repmat((ceil(self.Nx/2)+1)*ones(1,self.Ny),[1 self.nModes]),repmat(1:self.Ny,[1 self.nModes]),reshape(ones(1,self.Ny)'*(1:self.nModes),1,[]));
                nyquistIndicesForL = sub2ind(size(self.Omega),repmat(1:self.Nx,[1 self.nModes]),repmat((ceil(self.Ny/2)+1)*ones(1,self.Nx),[1 self.nModes]),reshape(ones(1,self.Nx)'*(1:self.nModes),1,[]));
                nyquistIndices = union(nyquistIndicesForK,nyquistIndicesForL);
            end
            
            % Sort the frequencies (for each mode) and distribute energy.
            % This algorithm is fairly complicated because we are using two
            % separate lists of frequencies: one for the gridded IW modes
            % and one for the 'external' modes.
            %
            % Note that it would appear that we are double counting the
            % number of waves at each frequency because we're included the
            % negative part of the hermitian conjugate. However, we also
            % have negative frequency waves, so this is justified.
            internalOmegaLinearIndices = reshape(1:numel(self.Omega),size(self.Omega));
            externalOmegaLinearIndices = 1:length(self.omega_ext);
            GM3Dint = zeros(size(self.Kh));
            GM3Dext = zeros(size(self.k_ext));
            indicesToSkip = reshape(internalOmegaLinearIndices(abs(self.Kh) < minK | abs(self.Kh) > maxK),[],1);
            for iMode = minMode:maxMode
                intOmegasLinearIndicesForIMode = reshape(internalOmegaLinearIndices(:,:,iMode),[],1); % Find the indices of everything with this iMode
                if excludeNyquist == 1
                    intOmegasLinearIndicesForIMode = setdiff( intOmegasLinearIndicesForIMode, nyquistIndices); % Now remove indices associated with the Nyquist
                end
                intOmegasLinearIndicesForIMode = setdiff( intOmegasLinearIndicesForIMode, indicesToSkip); % And remove indices outside the user specified wavenumber threshold.
                
                % Flatten the internal omegas (and their index)
                intOmegas = abs(self.Omega(intOmegasLinearIndicesForIMode));
                                
                % Now do the same for the external modes
                indices = find(self.j_ext == iMode);
                extOmegas = reshape(abs(self.omega_ext(indices)),[],1);
                extOmegasLinearIndicesForIMode = reshape(externalOmegaLinearIndices(indices),[],1);
                
                % Make a combined list, but note which list each omega came
                % from.
                if initializeModes == 0
                    allOmegas = cat(1,intOmegas,extOmegas);
                    allIndices = cat(1,intOmegasLinearIndicesForIMode,extOmegasLinearIndicesForIMode);
                    allSource = cat(1,zeros(size(intOmegas)), ones(size(extOmegas)));
                elseif initializeModes == 1
                    allOmegas = intOmegas;
                    allIndices = intOmegasLinearIndicesForIMode;
                    allSource = zeros(size(intOmegas));
                else
                    allOmegas = extOmegas;
                    allIndices = extOmegasLinearIndicesForIMode;
                    allSource = ones(size(extOmegas));
                end
                
                if isempty(allOmegas)
                    continue;
                end
                
                % Sort the frequencies for this mode.
                [sortedOmegas, sortedOmegasIndices] = sort(allOmegas);
                sortedIndices = allIndices(sortedOmegasIndices);
                sortedSource = allSource(sortedOmegasIndices);
                
                % Then find where the omegas differ.
                omegaDiffIndices = find(diff(sortedOmegas) > 0);               
                
                % Let's do a sanity check for users to make sure they don't
                % put too much energy in a single mode
                totalEnergyInThisMode = GM2D_int(self.f0,self.Nmax,iMode);
                
                lastIdx = 1;
                omega0 = sortedOmegas(lastIdx);
                leftDeltaOmega = 0;
                for idx=omegaDiffIndices'
                    currentIdx = idx+1;
                    nOmegas = currentIdx-lastIdx;
                                 
                    omega1 = sortedOmegas(idx + 1);
                    rightDeltaOmega = (omega1-omega0)/2;
                    
                    % This enforces our maximum allowed gap.
                    if rightDeltaOmega > maxDeltaOmega/2
                        rightDeltaOmega = maxDeltaOmega/2;
                    end
                    energyPerFrequency = GM2D_int(omega0-leftDeltaOmega,omega0+rightDeltaOmega,iMode)/nOmegas;
                    
                    if energyPerFrequency/totalEnergyInThisMode > energyWarningThreshold
                        warning('A j=%d mode has %d%% of the GM energy in a single mode',iMode,round(100*energyPerFrequency/totalEnergyInThisMode));
                    end
                    
                    for iIndex = lastIdx:(currentIdx-1)
                        if sortedSource(iIndex) == 0
                            GM3Dint(sortedIndices(iIndex)) = energyPerFrequency;
                        else
                            GM3Dext(sortedIndices(iIndex)) = energyPerFrequency;
                        end
                    end
                    
                    omega0 = omega1;
                    leftDeltaOmega = rightDeltaOmega;
                    lastIdx = currentIdx;
                end
                
                % Still have to deal with the last point.
                if lastIdx == 1
                    % There is only one point for this entire iMode!
                   leftDeltaOmega = omega0 - self.f0;
                   if leftDeltaOmega > maxDeltaOmega/2
                        leftDeltaOmega = maxDeltaOmega/2;
                   end
                   rightDeltaOmega = self.Nmax-omega0;
                   if rightDeltaOmega > maxDeltaOmega/2
                       rightDeltaOmega = maxDeltaOmega/2;
                   end
                else
                    % okay, so there's more than one point.
                    % so just be symmetric about this point,
                    rightDeltaOmega = leftDeltaOmega;
                    
                    % but don't let it exceed the buoyancy frequency
                    if omega0+rightDeltaOmega > self.Nmax
                        rightDeltaOmega = self.Nmax-omega0;
                    end
                    
                    % This enforces our maximum allowed gap.
                    if rightDeltaOmega > maxDeltaOmega/2
                        rightDeltaOmega = maxDeltaOmega/2;
                    end
                end
                
                nOmegas = length(sortedOmegas)+1-lastIdx;
                energyPerFrequency = GM2D_int(omega0-leftDeltaOmega,omega0+rightDeltaOmega,iMode)/nOmegas;
                if energyPerFrequency/totalEnergyInThisMode > energyWarningThreshold
                    warning('A j=%d mode has %d%% of the GM energy in a single mode',iMode,round(100*energyPerFrequency/totalEnergyInThisMode));
                end
                for iIndex = lastIdx:length(sortedOmegas)
                    if sortedSource(iIndex) == 0
                        GM3Dint(sortedIndices(iIndex)) = energyPerFrequency;
                    else
                        GM3Dext(sortedIndices(iIndex)) = energyPerFrequency;
                    end
                end
                
                
                
                
            end

            
            % At this stage GM3Dint contains all the energy, E_gm.
            % Now this needs to be split so that
            %       (1)     E<A_plus^2 + A_minus^2> = E_gm
            A = sqrt(GM3Dint/2); 
            
            % Each standard coefficient a(i,j,k) has equal conjugate, which
            % is already accounted for as having the same frequency. So
            % sum(a^2/2) really is the total energy.
            if shouldRandomizeAmplitude == 1
                A_plus = A.*GenerateHermitianRandomMatrix( size(self.K) );
                A_minus = A.*GenerateHermitianRandomMatrix( size(self.K) );
                
                self.U_ext = sqrt(2*GM3Dext./self.h_ext).*randn( size(self.h_ext) );
                self.PrecomputeExternalWaveCoefficients();                
            else
                % Randomize phases, but keep unit length
                A_plus = GenerateHermitianRandomMatrix( size(self.K) );
                A_minus = GenerateHermitianRandomMatrix( size(self.K) );
                
                goodIndices = abs(A_plus) > 0;
                A_plus(goodIndices) = A_plus(goodIndices)./abs(A_plus(goodIndices));
                A_plus = A.*A_plus;
                goodIndices = abs(A_minus) > 0;
                A_minus(goodIndices) = A_minus(goodIndices)./abs(A_minus(goodIndices));
                A_minus = A.*A_minus;
                
                % Check this factor of 2!!! Is the correct? squared
                % velocity to energy, I think.
                self.U_ext = sqrt(2*GM3Dext./self.h_ext);
                self.PrecomputeExternalWaveCoefficients();   
            end
                        
            A_minus(1,1,:) = conj(A_plus(1,1,:)); % Inertial motions go only one direction!
            
            self.GenerateWavePhases(A_plus,A_minus);
            
            self.didPreallocateAdvectionCoefficients = 0;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Add external waves to the model to fill out the spectrum
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = FillOutWaveSpectrum(self,maxTimeGap)
            % Add free/external waves to fill in the gaps of the gridded
            % solution. No gaps will be larger than 2*pi/maxTimeGap, and
            % the gaps will be smaller near f0.
            if nargin < 2
                maxTimeGap = 86140;
            end
            
            % the function dOmegas has an initial gap of dOmegaInitial, and
            % asymptotes to maxdOmega. Ramp up puts energy in the
            % lowest/most energetic modes.
            dOmegaInitial = 0.05*self.f0;
            maxdOmega = 2*pi/maxTimeGap; % Gaps will appear with observations longer than T = 2*pi/maxdOmega;
            Ln = -1/log(1-dOmegaInitial/maxdOmega);
            dOmegas = maxdOmega*(1-exp(-(1:100)'/Ln));
            gapOmegas = self.f0 + cumsum(dOmegas);
            omegaExt = [];
            jExt = [];
            for iMode = 1:self.nModes
                omegas = sort(reshape(abs(self.Omega(:,:,iMode)),[],1));
                
                % First fill in the lower triangle
                indices = find(gapOmegas < omegas(2));
                jExt = cat(1,jExt,iMode*ones(length(indices),1));
                omegaExt = cat(1,omegaExt,gapOmegas(indices));
                
                % Then fill in overly sized gaps
                diffOmega = diff(omegas);
                gapIndices = find(diffOmega>maxdOmega);
                for i=2:length(gapIndices)
                    n = ceil(diffOmega(gapIndices(i))/maxdOmega);
                    newOmegas = linspace(omegas(gapIndices(i)),omegas(gapIndices(i)+1),n+1)';
                    jExt = cat(1,jExt,iMode*ones(n-1,1));
                    omegaExt = cat(1,omegaExt,newOmegas(2:end-1));
                end
            end
            alphaExt = 2*pi*rand( size(omegaExt) );
            phiExt = 2*pi*rand( size(omegaExt) );
            UExt = zeros(size(omegaExt));
            
            self.SetExternalWavesWithFrequencies(omegaExt,alphaExt,jExt,phiExt,UExt,Normalization.kConstant);
            
            fprintf('Added %d external waves to fill out the GM spectrum.\n', length(omegaExt));
        end
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Eulerian---the dynamical fields on the grid at a given time
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [varargout] = VariableFieldsAtTime(self, t, varargin)
            % Primary method for accessing the dynamical variables on the
            % internal grid.
            %
            % Valid variable options are 'u', 'v', 'w', 'rho_prime', and
            % 'zeta'.
            varargout = cell(size(varargin));
            [varargout{:}] = self.InternalVariableFieldsAtTime(t, varargin{:});
            if ~isempty(self.k_ext)
                varargoutExt = cell(size(varargin));
                [varargoutExt{:}] = self.ExternalVariableFieldsAtTime(t, varargin{:});
                for iArg=1:length(varargout)
                    varargout{iArg} = varargout{iArg} + varargoutExt{iArg};
                end
            end
        end
        
        function [u,v,w] = VelocityFieldAtTime(self, t)
            % Return the velocity field, which is the sum of the gridded
            % and external/free waves at time t. Note that if you do not
            % need w, don't request it and it won't be computed.
            if nargout == 3
                [u,v,w] = self.VariableFieldsAtTime(t,'u','v','w');
            else
                [u,v] = self.VariableFieldsAtTime(t,'u','v');
            end
        end
        
        function rho = DensityFieldAtTime(self, t)
            % Return the density field, which is the sum of the density
            % mean field (variable in z) and the perturbation field
            % (variable in time and space).
            rho_bar = self.DensityMeanField;
            rho_prime = self.VariableFieldsAtTime(t,'rho_prime');
            rho = rho_bar + rho_prime;
        end
        
        function rho_bar = DensityMeanField(self)
            % Return the mean density field, which is a function of z only.
           rho_bar = self.RhoBarAtDepth(self.Z); 
        end
        
        function rho_prime = DensityPerturbationFieldAtTime(self, t)
            % Return the density perturbation field, which is computed as
            % the sum of gridded and external waves.
            rho_prime = self.VariableFieldsAtTime(t,'rho_prime');
        end
        
        function zeta = IsopycnalDisplacementFieldAtTime(self, t)
            % Returns the linear estimate of the isopycnal displacement,
            % commonly denote eta or zeta.
            zeta = self.VariableFieldsAtTime(t,'zeta');
        end
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Lagrangian---return the dynamical fields at a given location and time
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                       
        function [varargout] = VariablesAtTimePosition(self,t,x,y,z,interpolationMethod,varargin)
            % Primary method for accessing the dynamical variables on the
            % at any position or time.
            %
            % The method argument specifies how off-grid values should be
            % interpolated. Use 'exact' for the slow, but accurate,
            % spectral interpolation. Otherwise use 'spline' or some other
            % method used by Matlab's interp function.
            %
            % Valid variable options are 'u', 'v', 'w', 'rho_prime', and
            % 'zeta'.
            varargout = cell(size(varargin));
            [varargout{:}] = self.InternalVariablesAtTimePosition(t,x,y,z,interpolationMethod,varargin{:});
            if ~isempty(self.k_ext)
                varargoutExt = cell(size(varargin));
                [varargoutExt{:}] = self.ExternalVariablesAtTimePosition(t,x,y,z,varargin{:});
                for iArg=1:length(varargout)
                    varargout{iArg} = varargout{iArg} + varargoutExt{iArg};
                end
            end
        end
        
        function uvw = VelocityFieldAtTimePosition(self,t,xyz,interpolationMethod)
            % useful for integration methods where dy/dt is best given with
            % y as a single variable.
            [u,v,w] = self.VelocityAtTimePosition(t,xyz(:,1),xyz(:,2),xyz(:,3),interpolationMethod);
            uvw = cat(2,u,v,w);
        end
        
        function [u,v,w] = VelocityAtTimePosition(self,t,x,y,z,interpolationMethod)
            if nargout == 3
                [u,v,w] = self.VariablesAtTimePosition(t,x,y,z,interpolationMethod,'u','v','w');
            else
                [u,v] = self.VariablesAtTimePosition(t,x,y,z,interpolationMethod,'u','v');
            end
        end
                
        function rho = DensityAtTimePosition(self,t,x,y,z,interpolationMethod)
            rho_bar = self.RhoBarAtDepth(z);
            rho_prime = self.VariablesAtTimePosition(t,x,y,z,interpolationMethod,'rho_prime');
            rho = rho_bar + rho_prime;
        end
        
        function rho_bar = DensityMeanAtDepth(self, z)
            % Return the mean density field, which is a function of z only.
            rho_bar = self.RhoBarAtDepth(z);
        end
        
        function rho_prime = DensityPerturbationAtTimePosition(self, t,x,y,z,interpolationMethod)
            rho_prime = self.VariablesAtTimePosition(t,x,y,z,interpolationMethod,'rho_prime');
        end
        
        function zeta = IsopycnalDisplacementAtTimePosition(self,t,x,y,z,interpolationMethod)
            zeta = self.VariablesAtTimePosition(t,x,y,z,interpolationMethod,'zeta');
        end
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Internal/gridded wave modes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [varargout] = InternalVariableFieldsAtTime(self, t, varargin)
            % This is the primary function for computing the internal
            % gridded dynamical variables. It tries to be somewhat
            % optimized, by only computing the phase once, and only
            % computing the requested variables.
            %
            % Valid variable options are 'u', 'v', 'w', 'rho_prime', and
            % 'zeta'.
            if length(varargin) < 1
                return;
            end
            
            u = []; v = []; w = []; zeta = [];
            
            phase_plus = exp(sqrt(-1)*self.Omega*t);
            phase_minus = exp(-sqrt(-1)*self.Omega*t);
            
            varargout = cell(size(varargin));
            for iArg=1:length(varargin)
                if ( strcmp(varargin{iArg}, 'u') )
                    if isempty(u)
                        u_bar = self.u_plus.*phase_plus + self.u_minus.*phase_minus;
                        if self.performSanityChecks == 1
                            CheckHermitian(u_bar);
                        end
                        u = self.TransformToSpatialDomainWithF(u_bar);
                    end
                    varargout{iArg} = u;
                    
                    if ~isempty(self.B)
                       if isempty(self.u_g)
                          self.u_g = self.TransformToSpatialDomainWithF(-sqrt(-1)*(self.g/self.f0)*self.L.*self.B);
                       end
                       varargout{iArg} = varargout{iArg} + self.u_g;
                    end
                elseif ( strcmp(varargin{iArg}, 'v') )
                    if isempty(v)
                        v_bar = self.v_plus.*phase_plus + self.v_minus.*phase_minus;
                        if self.performSanityChecks == 1
                            CheckHermitian(v_bar);
                        end
                        v = self.TransformToSpatialDomainWithF(v_bar);
                    end
                    varargout{iArg} = v;
                    
                    if ~isempty(self.B)
                        if isempty(self.v_g)
                            self.v_g = self.TransformToSpatialDomainWithF(sqrt(-1)*(self.g/self.f0)*self.K.*self.B);
                        end
                        varargout{iArg} = varargout{iArg} + self.v_g;
                    end
                elseif ( strcmp(varargin{iArg}, 'w') )
                    if isempty(w)
                        w_bar = self.w_plus.*phase_plus + self.w_minus.*phase_minus;
                        if self.performSanityChecks == 1
                            CheckHermitian(w_bar);
                        end
                        w = self.TransformToSpatialDomainWithG(w_bar);
                    end
                    varargout{iArg} = w;
                elseif ( strcmp(varargin{iArg}, 'rho_prime') || strcmp(varargin{iArg}, 'zeta') )
                    if isempty(zeta)
                        zeta_bar = self.zeta_plus.*phase_plus + self.zeta_minus.*phase_minus;
                        if self.performSanityChecks == 1
                            CheckHermitian(zeta_bar);
                        end
                        zeta = self.TransformToSpatialDomainWithG(zeta_bar);
                        
                        if ~isempty(self.B)
                            if isempty(self.zeta_g)
                                self.zeta_g = self.TransformToSpatialDomainWithG(self.B);
                            end
                            zeta = zeta + self.zeta_g;
                        end
                    end
                    
                    if strcmp(varargin{iArg}, 'zeta')
                        varargout{iArg} = zeta;
                    elseif strcmp(varargin{iArg}, 'rho_prime')
                        varargout{iArg} = (self.rho0/9.81)*self.N2AtDepth(self.Z) .* zeta;
                    end
                else
                    error('Invalid option. You may request u, v, w, rho_prime, or zeta.');
                end
            end
        end
        
        function [varargout] = InternalVariablesAtTimePosition(self,t,x,y,z,method,varargin)
            % Returns the gridded/internal dynamical variables at any point in
            % space or time. Because these waves are gridded, we can use
            % interpolation (e.g, spline, linear) to deduce the variable
            % values off grid, or we can use spectral interpolation to get
            % 'exact' values.
            if self.advectionSanityCheck == 0
                self.advectionSanityCheck = 1;
                if (self.z(end)-self.z(1)) ~= self.Lz
                    warning('Vertical domain does not span the full depth of the ocean. This will lead to NaNs when advected particles leave the resolved domain.')
                end
            end
            
            varargout = cell(size(varargin));
            if strcmp(method,'exact')
                [varargout{:}] = self.InternalVariablesAtTimePositionExact(t,x,y,z,varargin{:});
            else
                [varargout{:}] = self.InternalVariableFieldsAtTime(t,varargin{:});
                [varargout{:}] = self.InterpolatedFieldAtPositionNewAndShiny(x,y,z,method,varargout{:});
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % External wave modes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [varargout] = ExternalVariableFieldsAtTime(self,t,varargin)
            % Returns the external wave modes at the grid points.
            varargout = cell(size(varargin));
            [varargout{:}] = self.ExternalVariablesAtTimePosition(t,reshape(self.X,[],1),reshape(self.Y,[],1), reshape(self.Z,[],1), varargin{:});
            for iArg=1:length(varargout)
                varargout{iArg} = reshape(varargout{iArg},self.Nx,self.Ny,self.Nz);
            end
        end
        
        function [varargout] = ExternalVariablesAtTimePosition(self,t,x,y,z,varargin)
            % This is the primary function for computing the external
            % dynamical variables. It tries to be somewhat
            % optimized, by only computing the phase once, and only
            % computing the requested variables.
            %
            % Valid variable options are 'u', 'v', 'w', 'rho_prime', and
            % 'zeta'.
            if (~isscalar(t) || ~iscolumn(x) || ~iscolumn(y) || ~iscolumn(z))
                error('t must be a scalar and (x,y,z) must be column vectors.\n');
            end
            isU = 0; isV = 0; isW = 0; isZeta= 0; isRho = 0;
            varargout = cell(size(varargin));
            for iArg=1:length(varargin)
                isU = isU | strcmp(varargin{iArg}, 'u');
                isV = isV | strcmp(varargin{iArg}, 'v');
                isW = isW | strcmp(varargin{iArg}, 'w');
                isZeta = isZeta | strcmp(varargin{iArg}, 'zeta');
                isRho = isRho | strcmp(varargin{iArg}, 'rho_prime');
            end
            
            if isU
                u = zeros(size(x));
            end
            if isV
                v=zeros(size(x));
            end
            if isW
                w=zeros(size(x));
            end
            if isZeta || isRho
                zeta = zeros(size(x));
            end
            
            for i=1:length(self.k_ext)
                % Compute the phase
                theta = x * self.k_ext(i) + y * self.l_ext(i) + (self.omega_ext(i)*t + self.phi_ext(i));
                
                % Compute the cosine & sine of the phase, if necessary
                if ( isU || isV || isZeta || isRho )
                    cos_theta = cos(theta);
                end
                if ( isU || isV || isW )
                    sin_theta = sin(theta);
                end
                
                % Compute the necessary vertical structure functions
                if ( isU || isV )
                    F = self.ExternalUVModeAtDepth(z,i);
                end
                if ( isW || isZeta || isRho )
                    G = self.ExternalWModeAtDepth(z,i);
                end
                
                % Now compute the requested variables
                if isU
                    u = u + (self.U_cos_ext(i) * cos_theta + self.U_sin_ext(i) * sin_theta).*F;
                end
                if isV
                    v = v + (self.V_cos_ext(i) * cos_theta + self.V_sin_ext(i) * sin_theta).*F;
                end
                if isW
                    w = w + (self.W_sin_ext(i) * sin_theta).*G;
                end
                if isZeta || isRho
                    zeta = zeta + (self.Zeta_cos_ext(i) * cos_theta) .* G;
                end
            end
            
            for iArg=1:length(varargin)
                if strcmp(varargin{iArg}, 'u')
                    varargout{iArg} = u;
                elseif strcmp(varargin{iArg}, 'v')
                    varargout{iArg} = v;
                elseif strcmp(varargin{iArg}, 'w')
                    varargout{iArg} = w;
                elseif strcmp(varargin{iArg}, 'zeta')
                    varargout{iArg} = zeta;
                elseif strcmp(varargin{iArg}, 'rho_prime')
                    varargout{iArg} = (self.rho0/self.g)*self.N2AtDepth(z) .* zeta;
                end
            end
        end
        
        function value = get.Rho_cos_ext(self)
            value = (self.rho0/self.g)*self.Zeta_cos_ext;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Useful methods for advecting particles
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [u] = VelocityAtTimePositionVector(self,t,p, interpolationMethod,shouldUseW)
            % Return the velocity at time t and position p, where size(p) =
            % [3 n]. The rows of p represent [x,y,z].
            % Optional argument is passed to the interpolation method. I
            % recommend either linear or spline.
            psize = size(p);
            if psize(1) == 3
                [u,v,w] = self.VelocityAtTimePosition(t,p(1,:),p(2,:),p(3,:), interpolationMethod);
                if exist('shouldUseW','var')
                    u = cat(1,u,v,shouldUseW.*w);
                else
                    u = cat(1,u,v,w);
                end
            else
                [u,v,w] = self.VelocityAtTimePosition(t,p(:,1),p(:,2),p(:,3), interpolationMethod);
                if exist('shouldUseW','var')
                    u = cat(2,u,v,shouldUseW.*w);
                else
                    u = cat(2,u,v,w);
                end
            end
        end
        
        function [u] = DrifterVelocityAtTimePositionVector(self,t,p, interpolationMethod)
            % Return the velocity at time t and position p, where size(p) =
            % [3 n]. The rows of p represent [x,y,z].
            psize = size(p);
            if psize(1) == 3
                [u,v,~] = self.VelocityAtTimePosition(t,p(1,:),p(2,:),p(3,:), interpolationMethod);
                u = cat(1,u,v,zeros(1,psize(2)));
            else
                [u,v,~] = self.VelocityAtTimePosition(t,p(:,1),p(:,2),p(:,3), interpolationMethod);
                u = cat(2,u,v,zeros(psize(1),1));
            end
        end
        
        function [zIsopycnal, rhoIsopycnal] = PlaceParticlesOnIsopycnal(self,x,y,z,varargin)
            % MAS 1/10/18 - added intext ('int' or 'both') to give option of using int vs. int+ext fields for rho_prime
            % Also added rhoIsopycnal as output.
            % Any floats with the same value of z will be moved to the same
            % isopycnal.
            %
            % interpolation should probably be 'spline'.
            % tolerance is in meters, 1e-8 should be fine.

            if mod(length(varargin),2) ~= 0
                error('Arguments must be given as name/value pairs.');
            end
            interpolationMethod = 'spline';
            tolerance = 1e-8;
            maxIterations = 200; % max number of iterations to attempt to converge
            useModes = 'all';
            shouldShowDiagnostics = 0;
            for iArg = 1:2:length(varargin)
                if strcmp(varargin{iArg}, 'interpolationMethod')
                    interpolationMethod = varargin{iArg+1};
                elseif strcmp(varargin{iArg}, 'tolerance')
                    tolerance = varargin{iArg+1};
                elseif strcmp(varargin{iArg}, 'maxIterations')
                    maxIterations = varargin{iArg+1};
                elseif strcmp(varargin{iArg}, 'shouldShowDiagnostics')
                    shouldShowDiagnostics = varargin{iArg+1};
                elseif strcmp(varargin{iArg}, 'useModes')
                    if strcmp(varargin{iArg+1}, 'all') || strcmp(varargin{iArg+1}, 'internalOnly') || strcmp(varargin{iArg+1}, 'externalOnly')
                        useModes = varargin{iArg+1};
                    else
                        error('Invalid option for initializeModes');
                    end
                else
                    error('Invalid argument');
                end
            end
            
            t = 0;            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % Wrap the particle positions, as necessary
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % (x,y) are periodic for the gridded solution
            x_tilde = mod(x,self.Lx);
            y_tilde = mod(y,self.Ly);
            
            dx = self.x(2)-self.x(1);
            x_index = floor(x_tilde/dx);
            dy = self.y(2)-self.y(1);
            y_index = floor(y_tilde/dy);
            
            % Identify the particles along the interpolation boundary
            if strcmp(interpolationMethod,'spline')
                S = 3+1; % cubic spline, plus buffer
            elseif strcmp(interpolationMethod,'linear')
                S = 1+1;
            end
            bp = x_index < S-1 | x_index > self.Nx-S | y_index < S-1 | y_index > self.Ny-S;
            
            % then do the same for particles along the boundary.
            x_tildeS = mod(x+S*dx,self.Lx);
            y_tildeS = mod(y+S*dy,self.Ly);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % Create the gridded internal density field rho and the interpolants
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Rho = self.InternalVariableFieldsAtTime(t,'rho_prime');
            RhoInterp = griddedInterpolant(self.X,self.Y,self.Z, Rho,interpolationMethod);
            if any(bp)
               RhoInterpS = griddedInterpolant(self.X,self.Y,self.Z, circshift(Rho,[S S 0]),interpolationMethod); 
            end
                    
            % Now let's place the floats along an isopycnal.
            zIsopycnal = z;
            zUnique = unique(z);
            
            if shouldShowDiagnostics == 1
                figure('Position',[100 100 700 600]);
            end
            
            % MAS 1/10/18: create scale factor (<=1.0) to avoid overshoot when trying to converge to isopycnal
            fac = 0.3;
            
            rho = zeros(size(zIsopycnal));                      % initialize rho array
            for zLevel = 1:length(zUnique)
                iterations = 0;
                zLevelIndices = (z==zUnique(zLevel));           % 1 if this is a zLevel, 0 if not
                
                nonboundaryIndices = zLevelIndices & ~bp;       % 1 if this is not a boundary index, 0 if it is
                % get rho for all non-boundary float locations
                rho(nonboundaryIndices) = RhoInterp(x_tilde(nonboundaryIndices),y_tilde(nonboundaryIndices),zIsopycnal(nonboundaryIndices));
                % now get rho for all boundary float locations
                if any(bp)
                    boundaryIndices = zLevelIndices & bp;
                    rho(boundaryIndices) = RhoInterpS(x_tildeS(boundaryIndices),y_tildeS(boundaryIndices),zIsopycnal(boundaryIndices));
                end
                
                % MAS 1/10/18: now make rho = gridded rho_prime + rho_bar (+ external rho_prime)
                if strcmp(useModes,'internalOnly')==1
                  rho(zLevelIndices) = rho(zLevelIndices) + self.RhoBarAtDepth(zIsopycnal(zLevelIndices));
                  elseif strcmp(useModes,'externalOnly')==1
                  rho(zLevelIndices) = self.RhoBarAtDepth(zIsopycnal(zLevelIndices)) + self.ExternalVariablesAtTimePosition(t,x(zLevelIndices),y(zLevelIndices),zIsopycnal(zLevelIndices),'rho_prime');
                elseif strcmp(useModes,'all')==1
                  rho(zLevelIndices) = rho(zLevelIndices) + self.RhoBarAtDepth(zIsopycnal(zLevelIndices)) + self.ExternalVariablesAtTimePosition(t,x(zLevelIndices),y(zLevelIndices),zIsopycnal(zLevelIndices),'rho_prime');
                end
                
                dRho = rho(zLevelIndices) - mean(rho(zLevelIndices));
                dz = dRho * 9.81./(self.N2AtDepth(zIsopycnal(zLevelIndices))*self.rho0);
                zIsopycnal(zLevelIndices) = zIsopycnal(zLevelIndices)+dz;
                
                while( max(abs(dz)) > tolerance && iterations < maxIterations )
                    rho(nonboundaryIndices) = RhoInterp(x_tilde(nonboundaryIndices),y_tilde(nonboundaryIndices),zIsopycnal(nonboundaryIndices));
                   if any(bp)
                        rho(boundaryIndices) = RhoInterpS(x_tildeS(boundaryIndices),y_tildeS(boundaryIndices),zIsopycnal(boundaryIndices));
                    end

                    % MAS 1/10/18: now make rho = gridded rho_prime + rho_bar (+ external rho_prime)
                    if strcmp(useModes,'internalOnly')==1
                        rho(zLevelIndices) = rho(zLevelIndices) + self.RhoBarAtDepth(zIsopycnal(zLevelIndices));
                    elseif strcmp(useModes,'externalOnly')==1
                        rho(zLevelIndices) = self.RhoBarAtDepth(zIsopycnal(zLevelIndices)) + self.ExternalVariablesAtTimePosition(t,x(zLevelIndices),y(zLevelIndices),zIsopycnal(zLevelIndices),'rho_prime');
                    elseif strcmp(useModes,'all')==1
                        rho(zLevelIndices) = rho(zLevelIndices) + self.RhoBarAtDepth(zIsopycnal(zLevelIndices)) + self.ExternalVariablesAtTimePosition(t,x(zLevelIndices),y(zLevelIndices),zIsopycnal(zLevelIndices),'rho_prime');
                    end

                    dRho = rho(zLevelIndices) - mean(rho(zLevelIndices));
                    dz = dRho * 9.81./(self.N2AtDepth(zIsopycnal(zLevelIndices))*self.rho0);
                    zIsopycnal(zLevelIndices) = zIsopycnal(zLevelIndices)+fac*dz;
                    iterations = iterations + 1;
                    if shouldShowDiagnostics == 1
                        fprintf('  PlaceParticlesOnIsopycnal: Iteration = %d. Mean dz=%3d m of isopycnal at z=%.1f m\n',iterations,mean(abs(dz)),mean(z(zLevelIndices)));
                        if iterations==1
                            subplot(2,1,1)
                            plot(zIsopycnal(zLevelIndices),'ro')
                            hold on
                            subplot(2,1,2)
                            plot(rho(zLevelIndices),'ro')
                            hold on
                        else
                            subplot(2,1,1)
                            plot(zIsopycnal(zLevelIndices),'bo')
                            subplot(2,1,2)
                            plot(rho(zLevelIndices),'bo')
                        end
                    end
                end
                % Do this one more time now that out of 'while' loop to get final rho - we'll pass this back along with zIsopycnal
                rho(nonboundaryIndices) = RhoInterp(x_tilde(nonboundaryIndices),y_tilde(nonboundaryIndices),zIsopycnal(nonboundaryIndices));
                if any(bp)
                  rho(boundaryIndices) = RhoInterpS(x_tildeS(boundaryIndices),y_tildeS(boundaryIndices),zIsopycnal(boundaryIndices));
                end

                % MAS 1/10/18: now make rho = gridded rho_prime + rho_bar (+ external rho_prime)
                if strcmp(useModes,'internalOnly')==1
                    rho(zLevelIndices) = rho(zLevelIndices) + self.RhoBarAtDepth(zIsopycnal(zLevelIndices));
                elseif strcmp(useModes,'externalOnly')==1
                    rho(zLevelIndices) = self.RhoBarAtDepth(zIsopycnal(zLevelIndices)) + self.ExternalVariablesAtTimePosition(t,x(zLevelIndices),y(zLevelIndices),zIsopycnal(zLevelIndices),'rho_prime');
                elseif strcmp(useModes,'all')==1
                  rho(zLevelIndices) = rho(zLevelIndices) + self.RhoBarAtDepth(zIsopycnal(zLevelIndices)) + self.ExternalVariablesAtTimePosition(t,x(zLevelIndices),y(zLevelIndices),zIsopycnal(zLevelIndices),'rho_prime');
                end

                if shouldShowDiagnostics == 1
                    subplot(2,1,1)
                    plot(zIsopycnal(zLevelIndices),'m*')
                    ylabel('float depth (m)')
                    subplot(2,1,2)
                    plot(rho(zLevelIndices),'m*')
                    xlabel('Float #')
                    ylabel('float \rho (kg/m^3)')
                end
                
                fprintf('Num iterations=%d. All floats within %.2g m of isopycnal at z=%.1f meters\n',iterations,max(abs(dz)),mean(z(zLevelIndices)) )
            end
            rhoIsopycnal = rho;
        end
        
        function ShowDiagnostics(self)
            % Display various diagnostics about the simulation.
            omega = abs(self.Omega);
            fprintf('Model resolution is %.2f x %.2f x %.2f meters.\n', self.x(2)-self.x(1), self.y(2)-self.y(1), self.z(2)-self.z(1));
            fprintf('The ratio Nmax/f0 is %.1f.\n', self.Nmax/self.f0);
            fprintf('Discretization effects will become apparent after %.1f hours in the frequency domain as the fastest modes traverse the domain.\n', max([self.Lx self.Ly])/max(max(max(self.C)))/3600);
            sortedOmega = sort(unique(reshape(omega(:,:,1),1,[])));
            fprintf('j=1 mode has discrete frequencies (%.4f f0, %.4f f0, ..., %.4f Nmax, %.4f Nmax)\n', sortedOmega(1)/self.f0, sortedOmega(2)/self.f0, sortedOmega(end-1)/self.Nmax, sortedOmega(end)/self.Nmax);
            dOmega = (sortedOmega(2)-sortedOmega(1))/2;
            T = 2*pi/dOmega;
            fprintf('The gap between these two lowest frequencies will be fully resolved after %.1f hours\n', T/3600);
            sortedOmega = sort(unique(reshape(omega(:,:,end),1,[])));
            fprintf('j=%d mode has discrete frequencies (%.4f f0, %.4f f0, ..., %.4f Nmax, %.4f Nmax)\n', self.nModes, sortedOmega(1)/self.f0, sortedOmega(2)/self.f0, sortedOmega(end-1)/self.Nmax, sortedOmega(end)/self.Nmax);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computes the phase information given the amplitudes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function GenerateWavePhases(self, U_plus, U_minus)
            % Given amplitudes A_plus,A_minus, this initializes the
            % amplitudes and phases for the dynamical variables. Generally
            % you don't want to call this directly, unless you have some
            % very specific use case.
            self.Amp_plus = U_plus; self.Amp_minus = U_minus;
            alpha = atan2(self.L,self.K);
            omega = abs(self.Omega); % The following definitions assume omega > 0.
            if abs(self.f0) < 1e-14 % This handles the f=0 case.
                omega(omega == 0) = 1;
            end
            fOverOmega = self.f0 ./ omega;
            
            % Without calling MakeHermitian, this doesn't deal with l=0.
            self.u_plus = U_plus .* InternalWaveModel.MakeHermitian( ( -sqrt(-1)*fOverOmega .* sin(alpha) + cos(alpha) )./sqrt(self.h) );
            self.u_minus = U_minus .* InternalWaveModel.MakeHermitian( (sqrt(-1)*fOverOmega .* sin(alpha) + cos(alpha) )./sqrt(self.h) );
            
            self.v_plus = U_plus .* InternalWaveModel.MakeHermitian( ( sqrt(-1)*fOverOmega .* cos(alpha) + sin(alpha) )./sqrt(self.h) );
            self.v_minus = U_minus .* InternalWaveModel.MakeHermitian( ( -sqrt(-1)*fOverOmega .* cos(alpha) + sin(alpha) )./sqrt(self.h) );
            
            self.w_plus = U_plus .* InternalWaveModel.MakeHermitian(-sqrt(-1) *  self.Kh .* sqrt(self.h) );
            self.w_minus = U_minus .* InternalWaveModel.MakeHermitian( -sqrt(-1) * self.Kh .* sqrt(self.h) );
            
            self.zeta_plus = U_plus .* InternalWaveModel.MakeHermitian( -self.Kh .* sqrt(self.h) ./ omega );
            self.zeta_minus = U_minus .* InternalWaveModel.MakeHermitian( self.Kh .* sqrt(self.h) ./ omega );
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Diagnostic shows a plot of all resolved wave modes, both gridded
        % and external
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        
        function omegaResolved = ResolvedFrequenciesAtMode(self,jMode)
            Om = self.Omega(:,:,jMode);
            OmNyquist = InternalWaveModel.NyquistWavenumbers(Om);  
            omegaResolved = sort(reshape(unique(abs(Om(~OmNyquist))),1,[]));
        end
        
        function ShowResolvedModesPlot(self,maxDeltaOmega)
            xAxisMax = 15*self.f0;
            omegaAxis = linspace(self.f0,xAxisMax,1000)';
            modeAxis = (1:15);
            if ~exist('maxDeltaOmega', 'var')
                maxDeltaOmega = self.Nmax-self.f0;
            end
            
            j_star = 3; % vertical mode roll-off
            L_gm = 1.3e3; % thermocline exponential scale, meters
            invT_gm = 5.2e-3; % reference buoyancy frequency, radians/seconds
            E_gm = 6.3e-5; % non-dimensional energy parameter
            E = L_gm*L_gm*L_gm*invT_gm*invT_gm*E_gm;
            N0 = invT_gm;
            
            H = (j_star+(1:3000)).^(-5/2);
            H_norm = 1/sum(H);
            B_norm = 1/atan(sqrt(N0*N0/(self.f0*self.f0)-1));
            
            GM2D_int = @(omega0,omega1,j) E*H_norm*B_norm*((j+j_star).^(-5/2))*(atan(self.f0/sqrt(omega0*omega0-self.f0*self.f0)) - atan(self.f0/sqrt(omega1*omega1-self.f0*self.f0)));
            GM2D_uv_int = @(omega0,omega1,j) E*H_norm*B_norm*((j+j_star).^(-5/2))*( self.f0*sqrt(omega1*omega1-self.f0*self.f0)/(2*omega1*omega1) - (3/2)*atan(self.f0/sqrt(omega1*omega1-self.f0*self.f0)) - self.f0*sqrt(omega0*omega0-self.f0*self.f0)/(2*omega0*omega0) + (3/2)*atan(self.f0/sqrt(omega0*omega0-self.f0*self.f0)));
            GM2D_w_int = @(omega0,omega1,j) E*H_norm*B_norm*((j+j_star).^(-5/2))*( self.f0*sqrt(omega1*omega1-self.f0*self.f0) + self.f0*self.f0*atan(self.f0/sqrt(omega1*omega1-self.f0*self.f0)) - self.f0*sqrt(omega0*omega0-self.f0*self.f0) - self.f0*self.f0*atan(self.f0/sqrt(omega0*omega0-self.f0*self.f0)));
            GM2D_zeta_int = @(omega0,omega1,j) E*H_norm*B_norm*((j+j_star).^(-5/2))*( ((omega1*omega1-self.f0*self.f0)^(3/2))/(2*self.f0*omega1*omega1) - (1/2)*atan(self.f0/sqrt(omega1*omega1-self.f0*self.f0)) - sqrt(omega1*omega1-self.f0*self.f0)/(2*self.f0) - ((omega0*omega0-self.f0*self.f0)^(3/2))/(2*self.f0*omega0*omega0) + (1/2)*atan(self.f0/sqrt(omega0*omega0-self.f0*self.f0)) + sqrt(omega0*omega0-self.f0*self.f0)/(2*self.f0) );
            
            TE = zeros(length(omegaAxis),length(modeAxis));
            HKE = zeros(length(omegaAxis),length(modeAxis));
            VKE = zeros(length(omegaAxis),length(modeAxis));
            IE = zeros(length(omegaAxis),length(modeAxis));
            for jMode = modeAxis
                for iOmega = 1:(length(omegaAxis)-1)
                    TE(iOmega,jMode) = GM2D_int(omegaAxis(iOmega),omegaAxis(iOmega+1),jMode);
                    HKE(iOmega,jMode) = GM2D_uv_int(omegaAxis(iOmega),omegaAxis(iOmega+1),jMode);
                    VKE(iOmega,jMode) = GM2D_w_int(omegaAxis(iOmega),omegaAxis(iOmega+1),jMode);
                    IE(iOmega,jMode) = GM2D_zeta_int(omegaAxis(iOmega),omegaAxis(iOmega+1),jMode);
                end
            end
            
            % This shift is applied to the points along the omega axis for visual aid.
            omega_epsilon = 0.0;
            
            ticks = linspace(self.f0,xAxisMax,5);
            
            labels = cell(length(ticks),1);
            labels{1} = 'f_0';
            for i=2:(length(ticks)-1)
                labels{i} = sprintf('%df_0',round(ticks(i)/self.f0));
            end
            if xAxisMax == N0
                labels{length(ticks)} = 'N_0';
            else
                labels{length(ticks)} = sprintf('%df_0',round(xAxisMax/self.f0));
            end
            
            vticks = (1.5:1:max(modeAxis))';
            vlabels = cell(length(vticks),1);
            for i=1:length(vticks)
                vlabels{i} = sprintf('%d',floor(vticks(i)));
            end
            
            scale = @(a) log10(a);
            
            figHandle = figure('Position',[100 100 700 600]);
            
            sp1 = subplot(2,2,1);
            pcolor(omegaAxis,modeAxis,scale(TE/max(max(TE)))'), hold on
            caxis([-2 0])
            shading flat
            for jMode = modeAxis
                omega = self.ResolvedFrequenciesAtMode(jMode);
                scatter(omega+omega_epsilon,jMode*ones(size(omega))+0.5,16*ones(size(omega)),'filled', 'MarkerFaceColor', 0*[1 1 1])
            end
            scatter(self.omega_ext+omega_epsilon,self.j_ext+0.5,16*ones(size(self.omega_ext)),'filled', 'MarkerFaceColor', 1*[1 1 1])
            for iMode = 1:self.nModes
                omega = abs(self.omega_ext(self.j_ext==iMode));
                omega = sort(unique(cat(2,omega,reshape(abs(self.Omega(:,:,iMode)),1,[]))));
                dOmega = diff(omega);
                indices = find(dOmega>maxDeltaOmega);
                for iBox=indices
                    rectangle('Position',[omega(iBox)+omega_epsilon+maxDeltaOmega/2 iMode dOmega(iBox)-maxDeltaOmega 1]);
                end
            end
            title('total')
            ylabel('vertical mode')
            xticks([])
            yticks(vticks)
            yticklabels(vlabels)
            
            sp2 = subplot(2,2,2);
            pcolor(omegaAxis,modeAxis,scale(HKE/max(max(HKE)))'), hold on
            caxis([-2 0])
            shading flat
            for jMode = modeAxis
                omega = self.ResolvedFrequenciesAtMode(jMode);
                scatter(omega+omega_epsilon,jMode*ones(size(omega))+0.5,16*ones(size(omega)),'filled', 'MarkerFaceColor', 0*[1 1 1])
            end
            scatter(self.omega_ext+omega_epsilon,self.j_ext+0.5,16*ones(size(self.omega_ext)),'filled', 'MarkerFaceColor', 1*[1 1 1])
            for iMode = 1:self.nModes
                omega = abs(self.omega_ext(self.j_ext==iMode));
                omega = sort(unique(cat(2,omega,reshape(abs(self.Omega(:,:,iMode)),1,[]))));
                dOmega = diff(omega);
                indices = find(dOmega>maxDeltaOmega);
                for iBox=indices
                    rectangle('Position',[omega(iBox)+omega_epsilon+maxDeltaOmega/2 iMode dOmega(iBox)-maxDeltaOmega 1]);
                end
            end
            title('horizontal')
            % ylabel('vertical mode')
            yticks([])
            % xlabel('frequency')
            xticks([])
            % xticklabels(labels)
            
            sp3 = subplot(2,2,3);
            pcolor(omegaAxis,modeAxis,scale(VKE/max(max(VKE)))'), hold on
            caxis([-2 0])
            shading flat
            for jMode = modeAxis
                omega = self.ResolvedFrequenciesAtMode(jMode);
                scatter(omega+omega_epsilon,jMode*ones(size(omega))+0.5,16*ones(size(omega)),'filled', 'MarkerFaceColor', 0*[1 1 1])
            end
            scatter(self.omega_ext+omega_epsilon,self.j_ext+0.5,16*ones(size(self.omega_ext)),'filled', 'MarkerFaceColor', 1*[1 1 1])
            for iMode = 1:self.nModes
                omega = abs(self.omega_ext(self.j_ext==iMode));
                omega = sort(unique(cat(2,omega,reshape(abs(self.Omega(:,:,iMode)),1,[]))));
                dOmega = diff(omega);
                indices = find(dOmega>maxDeltaOmega);
                for iBox=indices
                    rectangle('Position',[omega(iBox)+omega_epsilon+maxDeltaOmega/2 iMode dOmega(iBox)-maxDeltaOmega 1]);
                end
            end
            title('vertical')
            ylabel('vertical mode')
            xlabel('frequency')
            xticks(ticks)
            xticklabels(labels)
            yticks(vticks)
            yticklabels(vlabels)
            
            sp4 = subplot(2,2,4);
            pcolor(omegaAxis,modeAxis,scale(IE/max(max(IE)))'), hold on
            caxis([-2 0])
            shading flat
            for jMode = modeAxis
                omega = self.ResolvedFrequenciesAtMode(jMode);
                scatter(omega+omega_epsilon,jMode*ones(size(omega))+0.5,16*ones(size(omega)),'filled', 'MarkerFaceColor', 0*[1 1 1])
            end
            scatter(self.omega_ext+omega_epsilon,self.j_ext+0.5,16*ones(size(self.omega_ext)),'filled', 'MarkerFaceColor', 1*[1 1 1])
            for iMode = 1:self.nModes
                omega = abs(self.omega_ext(self.j_ext==iMode));
                omega = sort(unique(cat(2,omega,reshape(abs(self.Omega(:,:,iMode)),1,[]))));
                dOmega = diff(omega);
                indices = find(dOmega>maxDeltaOmega);
                for iBox=indices
                    rectangle('Position',[omega(iBox)+omega_epsilon+maxDeltaOmega/2 iMode dOmega(iBox)-maxDeltaOmega 1]);
                end
            end
            title('isopycnal')
            % ylabel('vertical mode')
            yticks([])
            xlabel('frequency')
            xticks(ticks)
            xticklabels(labels)
            % c = colorbar;
            
            dy = 0.045;
            sp1.Position(2) = sp1.Position(2) - dy;
            sp1.Position(4) = sp1.Position(4) + dy;
            sp3.Position(4) = sp3.Position(4) + dy;
            sp2.Position(2) = sp2.Position(2) - dy;
            sp2.Position(4) = sp2.Position(4) + dy;
            sp4.Position(4) = sp4.Position(4) + dy;
            
            dx = 0.045;
            sp1.Position(3) = sp1.Position(3) + dx;
            sp3.Position(3) = sp3.Position(3) + dx;
            sp2.Position(1) = sp2.Position(1) - dx;
            sp2.Position(3) = sp2.Position(3) + dx;
            sp4.Position(1) = sp4.Position(1) - dx;
            sp4.Position(3) = sp4.Position(3) + dx;
            
            figHandle.NextPlot = 'add';
            a = axes;
            
            %// Set the title and get the handle to it
            ht = title(sprintf('%dkm x %dkm x %dkm (%dx%dx%d)',round(self.Lx/1e3),round(self.Ly/1e3),round(self.Lz/1e3),self.Nx,self.Ny,self.Nz),'FontSize', 24);
            ht.Position = [0.5 1.04 0.5];
            
            %// Turn the visibility of the axes off
            a.Visible = 'off';
            
            %// Turn the visibility of the title on
            ht.Visible = 'on';
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Diagnostic shows a plot of all resolved wave modes, both gridded
        % and external
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function ShowResolvedWavenumbersPlot(self)              % MAS 12/15/17
            figHandle = figure('Position',[100 100 700 600]);

            m_ext = self.j_ext*pi/self.Lz;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % make some plots
            subplot(3,2,1)
            scatter(abs(self.Omega(:)/self.f0),self.J(:),'k.')
            hold on
            scatter(abs(self.omega_ext(:)/self.f0),self.j_ext(:),'r.')
            xlabel('\omega/f')
            ylabel('Vert. Mode')
            title('All Modes (int(k) + ext(r))')
            axis([1 max(abs(self.Omega(:)/self.f0)) 0.5+[0 self.nModes]])
            grid on

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            subplot(3,2,2)
            scatter3(self.Kh(:),self.M(:),abs(self.Omega(:)/self.f0),'k.')
            hold on
            scatter3(sqrt(self.k_ext(:).^2 + self.l_ext(:).^2),m_ext(:),abs(self.omega_ext(:)/self.f0),'rx')
            xlabel('K (rad/m)')
            ylabel('m (rad/m)')
            zlabel('\omega/f')
            title('All Modes (int(k) + ext(r))')
            axis tight

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            subplot(3,2,3)
            scatter(self.Kh(:),self.M(:),'k.')
            hold on
            scatter(sqrt(self.k_ext(:).^2 + self.l_ext(:).^2),m_ext(:),'r.')
            xlabel('K_h (rad/m)')
            ylabel('m (rad/m)')
            axis tight
            ax = axis;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            subplot(3,2,4)
            scatter(abs(self.K(:)),self.M(:),'k.')
            hold on
            scatter(abs(self.k_ext(:)),m_ext(:),'r.')
            xlabel('k (rad/m)')
            ylabel('m (rad/m)')
            axis tight
            ax = axis;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            subplot(3,2,5)
            scatter(abs(self.K(:)),abs(self.L(:)),'k.')
            hold on
            scatter(abs(self.k_ext(:)),abs(self.l_ext(:)),'r.')

            xlabel('k (rad/m)')
            ylabel('l (rad/m)')
            axis tight

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            subplot(3,2,6)
            scatter(abs(self.L(:)),self.M(:),'k.')
            hold on
            scatter(abs(self.l_ext(:)),m_ext(:),'r.')
            xlabel('l (rad/m)')
            ylabel('m (rad/m)')
            axis tight
            ax = axis;
        end

    end
    
    methods %(Access = protected)
        
        function self = GenerateWavePhasesForSpectralAdvection(self)
            % This mask zeros redundant hermitian conjugates, and doubles
            % the amplitude where appropriate, in order to extract the
            % correct amplitude for each wavenumber
            ConjugateMask = zeros(size(self.K));
            ConjugateMask(1:end,1:(self.Ny/2+1),1:end) = 2; % Primary conjugates contain half the amplitude
            ConjugateMask((self.Nx/2+1):end,1,:) = 0; % These guys are conjugate to (1:Nx/2,1,:)
            ConjugateMask(1,1,:) = 1; % self-conjugate
            ConjugateMask(self.Nx/2+1,1,:) = 1; % self-conjugate
            ConjugateMask(self.Nx/2+1,self.Ny/2+1,:) = 1; % self-conjugate
            ConjugateMask(1,self.Ny/2+1,:) = 1; % self-conjugate
            
            % First pull out the nonzero omega+ waves
            U_plus = (abs(self.Amp_plus)./sqrt(self.h)).*ConjugateMask;
            U_plus(1,1,:) = 2*U_plus(1,1,:);
            
            waveIndices = find( U_plus )';
            U = U_plus(waveIndices);
            self.h_int = self.h(waveIndices);
            self.phi_int = angle(self.Amp_plus(waveIndices));
            self.k_int = self.K(waveIndices);
            self.l_int = self.L(waveIndices);
            self.j_int = self.J(waveIndices);
            self.omega_int = self.Omega(waveIndices);
            
            % Then append the nonzero omega- waves
            U_minus = (abs(self.Amp_minus)./sqrt(self.h)).*ConjugateMask;
            U_minus(1,1,:) = 0*U_minus(1,1,:);
            
            waveIndices = find( U_minus )';
            U = cat(2, U, U_minus(waveIndices));
            self.h_int = cat(2,self.h_int,self.h(waveIndices));
            self.phi_int = cat(2, self.phi_int, angle(self.Amp_minus(waveIndices)));
            self.k_int = cat(2, self.k_int, self.K(waveIndices));
            self.l_int = cat(2, self.l_int, self.L(waveIndices));
            self.j_int = cat(2, self.j_int, self.J(waveIndices));
            self.omega_int = cat(2, self.omega_int, -self.Omega(waveIndices));
            
            alpha0 = atan2(self.l_int,self.k_int);
            Kh_ = sqrt( self.k_int .* self.k_int + self.l_int .* self.l_int);
            
            self.U_cos_int = U .* cos(alpha0);
            self.U_sin_int = U .* (self.f0 ./ self.omega_int) .* sin(alpha0);
            self.V_cos_int = U .* sin(alpha0);
            self.V_sin_int = -U .* (self.f0 ./ self.omega_int) .* cos(alpha0);
            self.W_sin_int = U .* Kh_ .* self.h_int;
            self.Zeta_cos_int = - U .* Kh_ .* self.h_int ./ self.omega_int;
            self.Rho_cos_int = - (self.rho0/9.81) .* U .* Kh_ .* self.h_int ./ self.omega_int;
            
            self.didPreallocateAdvectionCoefficients = 1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Begin initializing the wave field (internal use only)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = SetOmegaFromEigendepths(self, h)
            % Subclasses *must* call this method as part of intialization.
            self.h = h;
            self.C = sqrt( self.g*self.h );
            self.Omega = sqrt(self.C.*self.C.*self.K2 + self.f0*self.f0);         % Mode frequency
            
            % Create the hermitian conjugates of the phase vectors;
            if self.Ny > 1
                self.Omega(:,self.Ny/2+1:end,:) = -self.Omega(:,self.Ny/2+1:end,:);
            end
            if self.Nx > 1
                self.Omega((self.Nx/2+1):end,1,:) = -self.Omega((self.Nx/2+1):end,1,:);
            end
        end
               

                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % These functions return various associated with the external modes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % A few notes on speed. It's absolutely remarkable to me, but for
        % some reason, the following code is much slower than the
        % for-loops used in the actual algorithms.
        %
        % theta = x * self.k_ext + y * self.l_ext + (self.omega_ext*t + self.phi_ext); % [N M] + [1 M]
        % cos_theta = cos(theta); % [N M]
        % sin_theta = sin(theta); % [N M]
        % F = self.ExternalUVModesAtDepth(z); % [N M]
        % 
        % u = sum( (self.U_cos_ext .* cos_theta + self.U_sin_ext .* sin_theta) .* F, 2);
        % v = sum( (self.V_cos_ext .* cos_theta + self.V_sin_ext .* sin_theta) .* F, 2);
        %
        % This is easily tested with the following code:
        %         N = 1e6; % N points
        %         M = 1e3; % M waves
        %         
        %         x = rand(N,1);
        %         k = rand(1,M);
        %         A = rand(1,M);
        %         
        %         tic
        %         a = sum(A.*cos(x*k),2);
        %         toc
        %         
        %         tic
        %         a1 = cos(x*k)*A'; %  [N M] * [M 1]
        %         toc
        %         
        %         tic
        %         b = sum(A.*cos(bsxfun(@times,x,k)),2);
        %         toc
        %         
        %         c = zeros(size(x));
        %         tic
        %         for i=1:M
        %             c = c + A(i)*cos(k(i)*x);
        %         end
        %         toc
        
        function F = ExternalUVModeAtDepth(self, z, iWave)
            F = interp1(self.z,self.F_ext(:,iWave),z,'linear');
        end
        
        function G = ExternalWModeAtDepth(self, z, iWave)
            G = interp1(self.z,self.G_ext(:,iWave),z,'linear');
        end
        
   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % These functions return the exact/spectral velocity field for the
        % internal gridded field
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [varargout] = InternalVariablesAtTimePositionExact(self,t,x,y,z,varargin)
            % Return the velocity field associated with the gridded
            % velocity field, but using spectral interpolation, rather than
            % the FFT grid.
            % size(x) = [N 1]
            % size(phi) = [1 M]
            if self.didPreallocateAdvectionCoefficients == 0
                self.GenerateWavePhasesForSpectralAdvection();
            end
            
            isU = 0; isV = 0; isW = 0; isZeta= 0; isRho = 0;
            varargout = cell(size(varargin));
            for iArg=1:length(varargin)
                isU = isU | strcmp(varargin{iArg}, 'u');
                isV = isV | strcmp(varargin{iArg}, 'v');
                isW = isW | strcmp(varargin{iArg}, 'w');
                isZeta = isZeta | strcmp(varargin{iArg}, 'zeta');
                isRho = isRho | strcmp(varargin{iArg}, 'rho_prime');
            end
            
            if isU
                u = zeros(size(x));
            end
            if isV
                v=zeros(size(x));
            end
            if isW
                w=zeros(size(x));
            end
            if isZeta || isRho
                zeta = zeros(size(x));
            end
            
            for i=1:length(self.k_int)
                % Compute the phase
                theta = x * self.k_int(i) + y * self.l_int(i) + (self.omega_int(i)*t + self.phi_int(i));
                
                % Compute the cosine & sine of the phase, if necessary
                if ( isU || isV || isZeta || isRho )
                    cos_theta = cos(theta);
                end
                if ( isU || isV || isW )
                    sin_theta = sin(theta);
                end
                
                % Compute the necessary vertical structure functions
                if ( isU || isV )
                    F = self.InternalUVModeAtDepth(z,i);
                end
                if ( isW || isZeta || isRho )
                    G = self.InternalWModeAtDepth(z,i);
                end
                
                % Now compute the requested variables
                if isU
                    u = u + (self.U_cos_int(i) * cos_theta + self.U_sin_int(i) * sin_theta).*F;
                end
                if isV
                    v = v + (self.V_cos_int(i) * cos_theta + self.V_sin_int(i) * sin_theta).*F;
                end
                if isW
                    w = w + (self.W_sin_int(i) * sin_theta).*G;
                end
                if isZeta || isRho
                    zeta = zeta + (self.Zeta_cos_int(i) * cos_theta) .* G;
                end
            end
            
            for iArg=1:length(varargin)
                if strcmp(varargin{iArg}, 'u')
                    varargout{iArg} = u;
                elseif strcmp(varargin{iArg}, 'v')
                    varargout{iArg} = v;
                elseif strcmp(varargin{iArg}, 'w')
                    varargout{iArg} = w;
                elseif strcmp(varargin{iArg}, 'zeta')
                    varargout{iArg} = zeta;
                elseif strcmp(varargin{iArg}, 'rho_prime')
                    varargout{iArg} = (self.rho0/self.g)*self.N2AtDepth(z) .* zeta;
                end
            end
        end
                
        function varargout = InterpolatedFieldAtPosition(self,x,y,z,method,varargin)
            if nargin-5 ~= nargout
                error('You must have the same number of input variables as output variables');
            end
            
            % (x,y) are periodic for the gridded solution
            x_tilde = mod(x,self.Lx);
            y_tilde = mod(y,self.Ly);
            
            dx = self.x(2)-self.x(1);
            x_index = floor(x_tilde/dx);
            dy = self.y(2)-self.y(1);
            y_index = floor(y_tilde/dy);
            
            % Identify the particles along the interpolation boundary
            if strcmp(method,'spline')
                S = 3+1; % cubic spline, plus buffer
            elseif strcmp(method,'linear')
                S = 1+1;
            end
            bp = x_index < S-1 | x_index > self.Nx-S | y_index < S-1 | y_index > self.Ny-S;
            
            % then do the same for particles that along the boundary.
            x_tildeS = mod(x(bp)+S*dx,self.Lx);
            y_tildeS = mod(y(bp)+S*dy,self.Ly);
            
            varargout = cell(1,nargout);
            for i = 1:nargout
                U = varargin{i}; % gridded field
                u = zeros(size(x)); % interpolated value
                u(~bp) = interpn(self.X,self.Y,self.Z,U,x_tilde(~bp),y_tilde(~bp),z(~bp),method);
                if any(bp)
                    u(bp) = interpn(self.X,self.Y,self.Z,circshift(U,[S S 0]),x_tildeS,y_tildeS,z(bp),method);
                end
                varargout{i} = u;
            end     
        end
        
        function varargout = InterpolatedFieldAtPositionNewAndShiny(self,x,y,z,method,varargin)
            if nargin-5 ~= nargout
                error('You must have the same number of input variables as output variables');
            end
            
            % (x,y) are periodic for the gridded solution
            x_tilde = mod(x,self.Lx);
            y_tilde = mod(y,self.Ly);
            
            dx = self.x(2)-self.x(1);
            x_index = floor(x_tilde/dx);
            dy = self.y(2)-self.y(1);
            y_index = floor(y_tilde/dy);
            dz = self.z(2)-self.z(1);
            z_index = self.Nz+floor(z/dz)-1;
            
            if length(unique(diff(self.z)))>1
               error('This optimized interpolation does not yet work for uneven z') 
            end
                        
            % Identify the particles along the interpolation boundary
            if strcmp(method,'spline')
                S = 3+1; % cubic spline, plus buffer
            elseif strcmp(method,'linear')
                S = 1+1;
            end
            
            
            xrange = ((min(x_index)-S):(max(x_index)+S))+1;
            yrange = ((min(y_index)-S):(max(y_index)+S))+1;
            zrange = (max((min(z_index)-S),0):min((max(z_index)+S),self.Nz-1))+1;
                       
            varargout = cell(1,nargout);
            for i = 1:nargout
                U = varargin{i}; % gridded field
                u = interpn(self.X(xrange,yrange,zrange),self.Y(xrange,yrange,zrange),self.Z(xrange,yrange,zrange),U(xrange,yrange,zrange),x_tilde,y_tilde,z,method);
                varargout{i} = u;
            end
        end
                
    end
    
    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Forces a 3D matrix to be Hermitian, ready for an FFT (internal)
        %
        % The approach taken here is that the (k=-Nx/2..Nx/2,l=0..Ny/2+1) wave
        % numbers are primary, and the (k=-Nx/2..Nx/2,l=-Ny/2..1) are inferred as
        % conjugates. Also, the negative k wavenumbers for l=0. The Nyquist wave
        % numbers are set to zero to avoid complications.
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function A = MakeHermitian(A)
            M = size(A,1);
            N = size(A,2);
            K = size(A,3);
            
            % The order of the for-loop is chosen carefully here.
            for k=1:K
                for j=1:(N/2+1)
                    for i=1:M
                        ii = mod(M-i+1, M) + 1;
                        jj = mod(N-j+1, N) + 1;
                        if i == ii && j == jj
                            % A(i,j,k) = real(A(i,j,k)); % self-conjugate term
                            % This is not normally what you'd do, but we're being
                            % tricky by later adding the conjugate
                            if i == 1 && j == 1
                                continue;
                            else
                                A(i,j,k) = 0;
                            end
                        elseif j == N/2+1 % Kill the Nyquist, rather than fix it.
                            A(i,j,k) = 0;
                        else % we are letting l=0, k=Nx/2+1 terms set themselves again, but that's okay
                            A(ii,jj,k) = conj(A(i,j,k));
                        end
                    end
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Returns a matrix the same size as A with 1s at the 'redundant'
        % hermiation indices.
        function A = RedundantHermitianCoefficients(A)
            A = zeros(size(A));
            
            M = size(A,1);
            N = size(A,2);
            K = size(A,3);
            
            % The order of the for-loop is chosen carefully here.
            for k=1:K
                for j=1:(N/2+1)
                    for i=1:M
                        ii = mod(M-i+1, M) + 1;
                        jj = mod(N-j+1, N) + 1;
                        if i == ii && j == jj
                            if i == 1 && j == 1
                                continue;
                            else
                                A(i,j,k) = 0;
                            end
                        elseif j == N/2+1
                            A(i,j,k) = 0;
                        else
                            if j == 1 && i > M/2
                                continue;
                            end
                            A(ii,jj,k) = 1;
                        end
                    end
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Returns a matrix the same size as A with 1s at the Nyquist
        % frequencies.
        function A = NyquistWavenumbers(A)
            A = zeros(size(A));
            Nx = size(A,1);
            Ny = size(A,2);
            if Nx > 1
                A(ceil(Nx/2)+1,:) = 1;
            end
            if Ny > 1
                A(:,ceil(Ny/2)+1) = 1;
            end
        end
        
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Takes a Hermitian matrix and resets it back to the real amplitude.
function [A,phi,linearIndex] = ExtractNonzeroWaveProperties(Matrix)
M = size(Matrix,1);
N = size(Matrix,2);
K = size(Matrix,3);

A = [];
phi = [];
linearIndex = [];

% The order of the for-loop is chosen carefully here.
for k=1:K
    for j=1:(N/2+1)
        for i=1:M
            ii = mod(M-i+1, M) + 1;
            jj = mod(N-j+1, N) + 1;
            waveAmp = 0; wavePhase = 0;
            if i == ii && j == jj
                % self-conjugate term
                if i == 1 && j == 1
                    waveAmp = abs(Matrix(i,j,k));
                    wavePhase = angle(Matrix(i,j,k));
                else
                    continue;
                end
            elseif j == N/2+1 % Kill the Nyquist, rather than fix it.
                waveAmp = abs(Matrix(i,j,k));
                wavePhase = angle(Matrix(i,j,k));
            else % we are letting l=0, k=Nx/2+1 terms set themselves again, but that's okay 
%                 A(ii,jj,k) = conj(A(i,j,k));
                if j == 1 && i > M/2
                    continue;
                end
                waveAmp = 2*abs(Matrix(i,j,k));
                wavePhase = angle(Matrix(i,j,k));
            end
            if waveAmp > 0
               A = cat(1,A,waveAmp);
               phi = cat(1,phi,wavePhase);
               linearIndex = cat(1,linearIndex,sub2ind(size(Matrix),i,j,k));
            end
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Checks that the matrix is Hermitian.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = CheckHermitian(A)
M = size(A,1);
N = size(A,2);
K = size(A,3);

for k=1:K
   for i=M:-1:1
       for j=N:-1:1
           ii = mod(M-i+1, M) + 1;
           jj = mod(N-j+1, N) + 1;
           if A(i,j,k) ~= conj(A(ii,jj,k))
               error('Not hermitian conjugate')
           end
       end
   end
end

end

function A = GenerateHermitianRandomMatrix( size )

nX = size(1); nY = size(2); nZ = size(3);
A = InternalWaveModel.MakeHermitian(randn(size) + sqrt(-1)*randn(size) )/sqrt(2);
% A(1,1,:) = 2*real(A(1,1,:)); % Double the zero frequency
A(1,1,:) = 2*A(1,1,:); % Double the zero frequency
A(nX/2+1,1,:) = -2*real(A(nX/2+1,1,:)); % Double the Nyquist frequency
A(1,nY/2+1,:) = -2*real(A(1,nY/2+1,:)); % Double the Nyquist frequency
A(nX/2+1,nY/2+1,:) = -2*real(A(nX/2+1,nY/2+1,:)); % Double the Nyquist frequency
A(:,:,nZ) = zeros(nX,nY); % Because we can't resolve the last mode.

end


