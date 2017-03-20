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
    properties (Access = public)
        Lx, Ly, Lz % Domain size
        Nx, Ny, Nz % Number of points in each direction
        nModes
        latitude
        
        x, y, z
        k, l, j
        X,Y,Z
        K,L,J
        
        N2, Nmax
        
        K2, Kh, h, Omega, Omega_plus, Omega_minus, f0, C
        u_plus, u_minus, v_plus, v_minus, w_plus, w_minus, zeta_plus, zeta_minus
        
        Xh, Yh
        kExternal, lExternal, alphaExternal, omegaExternal, phiExternal, uExternal, FExternal, GExternal, hExternal
        
        version = 1.5
        performSanityChecks = 0
    end
    
    methods (Abstract, Access = protected)
        [F,G,h] = ModesAtWavenumber(self, k, norm ) % Return the normal modes and eigenvalue at a given wavenumber.
        [F,G,h] = ModesAtFrequency(self, omega, norm ) % Return the normal modes and eigenvalue at a given frequency.
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
            self.Nmax = max(N2);
            
            % Spectral domain, in radians
            dk = 1/self.Lx;          % fourier frequency
            self.k = 2*pi*([0:ceil(self.Nx/2)-1 -floor(self.Nx/2):-1]*dk)';
            
            dl = 1/self.Ly;          % fourier frequency
            self.l = 2*pi*([0:ceil(self.Ny/2)-1 -floor(self.Ny/2):-1]*dl)';
            
            self.j = (1:nModes)';
            
            [K,L,J] = ndgrid(self.k,self.l,self.j);
            [X,Y,Z] = ndgrid(self.x,self.y,self.z);
            [self.Xh,self.Yh] = ndgrid(self.x,self.y);
            
            self.L = L; self.K = K; self.J = J;
            self.X = X; self.Y = Y; self.Z = Z;
            
            self.f0 = 2 * 7.2921E-5 * sin( self.latitude*pi/180 );
            self.K2 = self.K.*self.K + self.L.*self.L;   % Square of the horizontal wavenumber
            self.Kh = sqrt(self.K2);
            
            nothing = zeros(size(self.K));
            self.h = nothing; self.Omega = nothing; self.Omega_plus = nothing;
            self.Omega_minus = nothing; self.C = nothing; self.u_plus = nothing;
            self.u_minus = nothing; self.v_plus = nothing; self.v_minus = nothing;
            self.w_plus = nothing; self.w_minus = nothing; self.zeta_plus = nothing;
            self.zeta_minus = nothing;
        end
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Create a single wave (public)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function period = InitializeWithPlaneWave(self, k0, l0, j0, UAmp, sign)
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
                sign=1;
            elseif l0 == 0 && k0 < 0
                k0 = -k0;
                sign = -1*sign;
                UAmp = -1*UAmp;
            elseif l0 < 0
                l0 = -l0;
                k0 = -k0;
                sign = -1*sign;
                UAmp = -1*UAmp;
            end
            
            % Rewrap (k0,l0) to follow standard FFT wrapping. l0 should
            % already be correct.
            if (k0 < 0)
                k0 = self.Nx + k0;
            end
                            
            ratio = self.UmaxGNormRatioForWave(k0, l0, j0);
            
            U = zeros(size(self.K));
            U(k0+1,l0+1,j0) = UAmp*ratio/2;
            if sign > 0
                A_plus = MakeHermitian(U);
                A_minus = zeros(size(U));
                A_minus(1,1,:) = A_plus(1,1,:); % Inertial oscillations are created using this trick.                
            else
                A_plus = zeros(size(U));
                A_minus = MakeHermitian(U);
            end
            
            self.GenerateWavePhases(A_plus,A_minus);
            
            period = 2*pi/self.Omega(k0+1,l0+1,j0);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Create a full Garrett-Munk spectrum (public)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function InitializeWithGMSpectrum(self, Amp)
            % GM Parameters
            j_star = 3;
            L_gm = 1.3e3; % thermocline exponential scale, meters
            invT_gm = 5.2e-3; % reference buoyancy frequency, radians/seconds
            E_gm = 6.3e-5; % non-dimensional energy parameter
            E = L_gm*L_gm*L_gm*invT_gm*invT_gm*E_gm*Amp;
            E = E*(self.Lz/L_gm); % This correction fixes the amplitude so that the HKE variance at a given depth matches (instead of depth integrated energy)
                      
            % Compute the proper vertical function normalization
            H = (j_star+(1:1024)).^(-5/2);
            H_norm = 1/sum(H);
            
            % Do the same for the frequency function.
            B_norm = 1/atan(sqrt(self.N0*self.N0/(self.f0*self.f0)-1));
            
            % This function tells you how much energy you need between two
            % frequencies for a given vertical mode.
            GM2D_int = @(omega0,omega1,j) E*H_norm*B_norm*((j+j_star).^(-5/2))*(atan(self.f0/sqrt(omega0*omega0-self.f0*self.f0)) - atan(self.f0/sqrt(omega1*omega1-self.f0*self.f0)));
            
            GM2D_uv_int = @(omega0,omega1,j) E*H_norm*B_norm*((j+j_star).^(-5/2))*( self.f0*sqrt(omega1*omega1-self.f0*self.f0)/(2*omega1*omega1) - (3/2)*atan(self.f0/sqrt(omega1*omega1-self.f0*self.f0)) - self.f0*sqrt(omega0*omega0-self.f0*self.f0)/(2*omega0*omega0) + (3/2)*atan(self.f0/sqrt(omega0*omega0-self.f0*self.f0)));
            GM2D_w_int = @(omega0,omega1,j) E*H_norm*B_norm*((j+j_star).^(-5/2))*( self.f0*sqrt(omega1*omega1-self.f0*self.f0) + self.f0*self.f0*atan(self.f0/sqrt(omega1*omega1-self.f0*self.f0)) - self.f0*sqrt(omega0*omega0-self.f0*self.f0) - self.f0*self.f0*atan(self.f0/sqrt(omega0*omega0-self.f0*self.f0)));
            GM2D_zeta_int = @(omega0,omega1,j) E*H_norm*B_norm*((j+j_star).^(-5/2))*( ((omega1*omega1-self.f0*self.f0)^(3/2))/(2*self.f0*omega1*omega1) - (1/2)*atan(self.f0/sqrt(omega1*omega1-self.f0*self.f0)) - sqrt(omega1*omega1-self.f0*self.f0)/(2*self.f0) - ((omega0*omega0-self.f0*self.f0)^(3/2))/(2*self.f0*omega0*omega0) + (1/2)*atan(self.f0/sqrt(omega0*omega0-self.f0*self.f0)) + sqrt(omega0*omega0-self.f0*self.f0)/(2*self.f0) );
            
            % Do a quick check to see how much energy is lost due to
            % limited vertical resolution.
            totalEnergy = 0;
            for mode=1:(max(self.j)/1)
                totalEnergy = totalEnergy + GM2D_int(self.f0,self.N0,mode);
            end
            fprintf('You are missing %.2f%% of the energy due to limited vertical modes.\n',100-100*totalEnergy/E);
            
            % Find the *second* lowest frequency
            [sortedOmegas, indices] = sort(reshape(abs(self.Omega(:,:,max(self.j)/2)),1,[]));
%             omegaStar = sortedOmegas(2);
            omegaStar = 1.6*self.f0;
            
            wVariancePerMode = [];
            for mode=1:(max(self.j)/2)
                wVariancePerMode(mode) = GM2D_w_int(self.f0+(min(min(self.Omega(2:end,2:end,mode)))-self.f0)/2,max(max(self.Omega(:,:,mode))),1);
                wVariancePerModeStar(mode) = GM2D_w_int(self.f0+(omegaStar-self.f0)/2,max(max(self.Omega(:,:,mode))),1);
            end
            
            shouldUseOmegaStar = 1;
            shouldUseMaxOmega = 0;
            
            % Sort the frequencies (for each mode) and distribute energy.
            GM3D = zeros(size(self.Kh));
            for iMode = 1:(max(self.j)/1)
                % Stride to the linear index for the full 3D matrix
                modeStride = (iMode-1)*size(self.Omega,1)*size(self.Omega,2);
                
                % Sort the linearized frequencies for this mode.
                [sortedOmegas, indices] = sort(reshape(abs(self.Omega(:,:,iMode)),1,[]));
                
                % Then find where the omegas differ.
                omegaDiffIndices = find(diff(sortedOmegas) > 0);
            
                lastIdx = 1;
                omega0 = sortedOmegas(lastIdx);
                leftDeltaOmega = 0;
                for idx=omegaDiffIndices
                    currentIdx = idx+1;
                    nOmegas = currentIdx-lastIdx;
                    
%                     if omega0 ~= obj.f0
%                         continue;
%                     end
                    
                    if shouldUseOmegaStar && omega0 == self.f0
                        omega1 = min(omegaStar,sortedOmegas(idx + 1));
                    elseif shouldUseMaxOmega == 1
                        omega1 = min(2.0*omega0,sortedOmegas(idx + 1));
                    else
                        omega1 = sortedOmegas(idx + 1);
                    end
                    rightDeltaOmega = (omega1-omega0)/2;
                      energyPerFrequency = GM2D_int(omega0-leftDeltaOmega,omega0+rightDeltaOmega,iMode)/nOmegas;
%                     energyPerFrequency = (GM2D_uv_int(omega0-leftDeltaOmega,omega0+rightDeltaOmega,iMode)/nOmegas)*(omega0*omega0/(omega0*omega0 + obj.f0*obj.f0));
                    
%                     if omega0 == obj.f0
%                         energyPerFrequency = GM2D_int(omega0-leftDeltaOmega,omega0+rightDeltaOmega,iMode)/nOmegas;
%                     else 
%                         energyPerFrequency = (GM2D_zeta_int(omega0-leftDeltaOmega,omega0+rightDeltaOmega,iMode)/nOmegas)*(omega0*omega0/(omega0*omega0 - obj.f0*obj.f0));
%                         energyPerFrequency = (GM2D_w_int(omega0-leftDeltaOmega,omega0+rightDeltaOmega,iMode)/nOmegas)*(1/(omega0*omega0 - obj.f0*obj.f0));
%                     end
                    
%                     energyPerFrequency = (GM2D_uv_int(omega0-leftDeltaOmega,omega0+rightDeltaOmega,iMode)/nOmegas)*(omega0*omega0/(omega0*omega0 + obj.f0*obj.f0));
                    
                    GM3D(indices(lastIdx:(currentIdx-1))+modeStride) = energyPerFrequency;
                    
                    if shouldUseOmegaStar && omega0 == self.f0
                        leftDeltaOmega = sortedOmegas(idx + 1) - (omega0+rightDeltaOmega);
                        omega0 = sortedOmegas(idx + 1);
                    elseif shouldUseMaxOmega == 1
                        omega0 = sortedOmegas(idx + 1);
                        leftDeltaOmega = rightDeltaOmega;
                    else
                        omega0 = omega1;
                        leftDeltaOmega = rightDeltaOmega;
                    end
                    
%                     omega0 = sortedOmegas(idx + 1); % same as omega0 = omega1, except when omega0 == f0
                     lastIdx = currentIdx;
%                     leftDeltaOmega = rightDeltaOmega;
                end
                % Still have to deal with the last point.
            end
            fprintf('After distributing energy across frequency and mode, you still have %.2f%% of reference GM energy.\n',100*sum(sum(sum(GM3D)))/E);
            fprintf('Due to restricted domain size, the j=1,k=l=0 mode contains %.2f%% the total energy.\n',100*GM3D(1,1,1)/sum(sum(sum(GM3D))));
            
            A = sqrt(GM3D/2); % Now split this into even and odd.
            
            % Randomize phases, but keep unit length
            A_plus = GenerateHermitianRandomMatrix( size(self.K) );
            A_minus = GenerateHermitianRandomMatrix( size(self.K) );
            
            goodIndices = abs(A_plus) > 0;
            A_plus(goodIndices) = A_plus(goodIndices)./abs(A_plus(goodIndices));
            A_plus = A.*A_plus;
            goodIndices = abs(A_minus) > 0;
            A_minus(goodIndices) = A_minus(goodIndices)./abs(A_minus(goodIndices));
            A_minus = A.*A_minus;
            
%              A_plus = A;
%              A_minus = A;        
            
%             A_plus = A.*GenerateHermitianRandomMatrix( size(obj.K) );
%             A_minus = A.*GenerateHermitianRandomMatrix( size(obj.K) );
            A_minus(1,1,:) = conj(A_plus(1,1,:)); % Intertial motions go only one direction!
            
            GM_sum = sum(sum(sum(GM3D)))/E;
            GM_random_sum = sum(sum(sum(A_plus.*conj(A_plus) + A_minus.*conj(A_minus)  )))/E;
            fprintf('The coefficients sum to %.2f%% GM given the scales, and the randomized field sums to %.2f%% GM\n', 100*GM_sum, 100*GM_random_sum);
            
            self.GenerateWavePhases(A_plus,A_minus);
        end
               

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Compute the dynamical fields at a given time (public)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [u,v] = VelocityFieldAtTime(self, t)
            phase_plus = exp(sqrt(-1)*self.Omega*t);
            phase_minus = exp(-sqrt(-1)*self.Omega*t);
            u_bar = self.u_plus.*phase_plus + self.u_minus.*phase_minus;
            v_bar = self.v_plus.*phase_plus + self.v_minus.*phase_minus;
            
            if self.performSanityChecks == 1
                CheckHermitian(u_bar);CheckHermitian(v_bar);
            end
            
            u = self.TransformToSpatialDomainWithF(u_bar);
            v = self.TransformToSpatialDomainWithF(v_bar);
            
            % Add the external waves to the solution
            if ~isempty(self.uExternal)
               [u_ext, v_ext] = self.ExternalVelocityFieldsAtTime(t);
               u = u + u_ext;
               v = v + v_ext;
            end
        end
        
        function [w,zeta] = VerticalFieldsAtTime(self, t)
            phase_plus = exp(sqrt(-1)*self.Omega*t);
            phase_minus = exp(-sqrt(-1)*self.Omega*t);
            w_bar = self.w_plus.*phase_plus + self.w_minus.*phase_minus;
            zeta_bar = self.zeta_plus.*phase_plus + self.zeta_minus.*phase_minus;
            
            if self.performSanityChecks == 1
                CheckHermitian(w_bar);CheckHermitian(zeta_bar);
            end
            
            w = self.TransformToSpatialDomainWithG(w_bar);
            zeta = self.TransformToSpatialDomainWithG(zeta_bar);
            
            if ~isempty(self.uExternal)
                [w_ext, zeta_ext] = self.ExternalVerticalFieldsAtTime(t);
                w = w + w_ext;
                zeta = zeta + zeta_ext;
            end
        end
        
        function omega = SetExternalWavesWithWavenumbers(self, k, l, j, phi, A, type)
            self.kExternal = k;
            self.lExternal = l;
            self.alphaExternal = atan2(l,k);
            self.phiExternal = phi;
            self.uExternal = A;
            
            if strcmp(type, 'energyDensity')
                norm = 'const_G_norm';
            elseif strcmp(type, 'maxU')
                norm = 'max_u';
            end
            
            g = 9.81;
            for iWave=1:length(j)
                K2h = k(iWave)*k(iWave) + l(iWave)*l(iWave);
                [F_ext,G_ext,h_ext] = self.ModesAtWavenumber(sqrt(K2h), norm);
                self.FExternal(:,iWave) = F_ext(:,j);
                self.GExternal(:,iWave) = G_ext(:,j);
                self.hExternal(iWave) = h_ext(j);
                self.omegaExternal(iWave) = sqrt( g*h_ext(iWave) * K2h + self.f0*self.f0 );
                if strcmp(type, 'energyDensity')
                    self.uExternal(iWave) = self.uExternal(iWave)/sqrt(h_ext(iWave));
                end
            end
            
            omega = self.omegaExternal;
        end
        
        function k = SetExternalWavesWithFrequencies(self, omega, alpha, j, phi, A, type)
            self.omegaExternal = omega;
            self.alphaExternal = alpha;
            self.phiExternal = phi;
            self.uExternal = A;
            
            if strcmp(type, 'energyDensity')
                norm = 'const_G_norm';
            elseif strcmp(type, 'maxU')
                norm = 'max_u';
            end
            
            g = 9.81;
            for iWave=1:length(j)
                [F_ext,G_ext,h_ext] = self.ModesAtFrequency(omega, norm);
                self.FExternal(:,iWave) = F_ext(:,j);
                self.GExternal(:,iWave) = G_ext(:,j);
                self.hExternal(iWave) = h_ext(j);
                K_horizontal = sqrt((omega*omega - self.f0*self.f0)/(g*h_ext(j)));
                self.kExternal(iWave) = K_horizontal*cos(alpha(iWave));
                self.lExternal(iWave) = K_horizontal*sin(alpha(iWave));
                if strcmp(type, 'energyDensity')
                    self.uExternal(iWave) = self.uExternal(iWave)/sqrt(h_ext(iWave));
                end
            end
            
            k = sqrt(self.kExternal*self.kExternal + self.lExternal*self.lExternal);
        end
    end
    
    methods (Access = protected)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Begin initializing the wave field (internal use only)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = SetOmegaFromEigendepths(self, h)
            % Subclasses *must* call this method as part of intialization.
            g = 9.81;
            self.h = h;
            self.C = sqrt( g*self.h );
            self.Omega = sqrt(self.C.*self.C.*self.K2 + self.f0*self.f0);         % Mode frequency
            
            % Create the hermitian conjugates of the phase vectors;
            self.Omega(:,(self.Ny/2+1):end,:) = -self.Omega(:,(self.Ny/2+1):end,:);
            self.Omega((self.Nx/2+1):end,1,:) = -self.Omega((self.Nx/2+1):end,1,:);
        end
        
        function ShowDiagnostics(self)
            omega = abs(self.Omega);
            fprintf('Model resolution is %.2f x %.2f x %.2f meters.\n', self.x(2)-self.x(1), self.y(2)-self.y(1), self.z(2)-self.z(1));
            fprintf('The ratio Nmax/f0 is %.1f.\n', self.Nmax/self.f0);
            fprintf('Discretization effects will become apparent after %.1f hours in the frequency domain as the fastest modes traverse the domain.\n', max([self.Lx self.Ly])/max(max(max(self.C)))/3600);
            sortedOmega = sort(unique(reshape(omega(:,:,1),1,[])));
            fprintf('j=1 mode has discrete frequencies (%.4f f0, %.4f f0, ..., %.4f N0, %.4f N0)\n', sortedOmega(1)/self.f0, sortedOmega(2)/self.f0, sortedOmega(end-1)/self.Nmax, sortedOmega(end)/self.Nmax);
            dOmega = (sortedOmega(2)-sortedOmega(1))/2;
            T = 2*pi/dOmega;
            fprintf('The gap between these two lowest frequencies will be fully resolved after %.1f hours\n', T/3600);
            sortedOmega = sort(unique(reshape(omega(:,:,end),1,[])));
            fprintf('j=%d mode has discrete frequencies (%.4f f0, %.4f f0, ..., %.4f N0, %.4f N0)\n', self.nModes, sortedOmega(1)/self.f0, sortedOmega(2)/self.f0, sortedOmega(end-1)/self.Nmax, sortedOmega(end)/self.Nmax);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computes the phase information given the amplitudes (internal)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function GenerateWavePhases(self, U_plus, U_minus)
            alpha = atan2(self.L,self.K);
            omega = abs(self.Omega); % The following definitions assume omega > 0.
            denominator = omega.*sqrt(self.h);
            
            % Without calling MakeHermitian, this doesn't deal with l=0.
            self.u_plus = U_plus .* MakeHermitian( ( -sqrt(-1)*self.f0 .* sin(alpha) + omega .* cos(alpha) )./denominator );
            self.u_minus = U_minus .* MakeHermitian( (sqrt(-1)*self.f0 .* sin(alpha) + omega .* cos(alpha) )./denominator );
            
            self.v_plus = U_plus .* MakeHermitian( ( sqrt(-1)*self.f0 .* cos(alpha) + omega .* sin(alpha) )./denominator );
            self.v_minus = U_minus .* MakeHermitian( ( -sqrt(-1)*self.f0 .* cos(alpha) + omega .* sin(alpha) )./denominator );
            
            self.w_plus = U_plus .* MakeHermitian(-sqrt(-1) *  self.Kh .* sqrt(self.h) );
            self.w_minus = U_minus .* MakeHermitian( -sqrt(-1) * self.Kh .* sqrt(self.h) );
            
            self.zeta_plus = U_plus .* MakeHermitian( -self.Kh .* sqrt(self.h) ./ omega );
            self.zeta_minus = U_minus .* MakeHermitian( self.Kh .* sqrt(self.h) ./ omega );
        end
        
        function [u,v] = ExternalVelocityFieldsAtTime(self, t)
            % Return the velocity field associated with the manually added
            % waves.
            nExternalWaves = length(self.uExternal);
            
            u = zeros(size(self.X));
            v = zeros(size(self.X));
            for iWave=1:nExternalWaves
                % Compute the two-dimensional phase vector (note we're
                % using Xh,Yh---the two-dimensional versions.
                k0 = self.kExternal(iWave);
                l0 = self.lExternal(iWave);
                omega0 = self.omegaExternal(iWave);
                phi0 = self.phiExternal(iWave);
                alpha0 = self.alphaExternal(iWave);
                U = self.uExternal(iWave);
                F = permute(self.FExternal(:,iWave),[3 2 1]);
                
                theta = k0 * self.Xh + l0 * self.Yh + omega0*t + phi0;
                
                % This .* should take [Nx Ny 1] and multiply by [1 1 Nz]
                u_wave = U*( cos(alpha0)*cos(theta) + (self.f0/omega0)*sin(alpha0)*sin(theta) ) .* F;
                v_wave = U*( sin(alpha0)*cos(theta) - (self.f0/omega0)*cos(alpha0)*sin(theta) ) .* F;
                
                u = u + u_wave;
                v = v + v_wave;
            end
        end
        
        function [w,zeta] = ExternalVerticalFieldsAtTime(self, t)
            % Return the velocity field associated with the manually added
            % waves.
            nExternalWaves = length(self.uExternal);
            
            w = zeros(size(self.X));
            zeta = zeros(size(self.X));
            for iWave=1:nExternalWaves
                % Compute the two-dimensional phase vector (note we're
                % using Xh,Yh---the two-dimensional versions.
                k0 = self.kExternal(iWave);
                l0 = self.lExternal(iWave);
                omega0 = self.omegaExternal(iWave);
                phi0 = self.phiExternal(iWave);
                h0 = self.hExternal(iWave);
                U = self.uExternal(iWave);
                G = permute(self.GExternal(:,iWave),[3 2 1]);
                
                theta = k0 * self.Xh + l0 * self.Yh + omega0*t + phi0;
                
                % This .* should take [Nx Ny 1] and multiply by [1 1 Nz]
                w_wave = U * sqrt(k0*k0+l0*l0) * h0 * sin(theta) .* G;
                zeta_wave = -U * sqrt(k0*k0+l0*l0) * (h0/omega0) * cos(theta) .* G;
                
                w = w + w_wave;
                zeta = zeta + zeta_wave;
            end
        end
    end
end

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
                % This is normally what you'd do, but we're being tricky
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
A = MakeHermitian(randn(size) + sqrt(-1)*randn(size) )/sqrt(2);
% A(1,1,:) = 2*real(A(1,1,:)); % Double the zero frequency
A(1,1,:) = 2*A(1,1,:); % Double the zero frequency
A(nX/2+1,1,:) = -2*real(A(nX/2+1,1,:)); % Double the Nyquist frequency
A(1,nY/2+1,:) = -2*real(A(1,nY/2+1,:)); % Double the Nyquist frequency
A(nX/2+1,nY/2+1,:) = -2*real(A(nX/2+1,nY/2+1,:)); % Double the Nyquist frequency
A(:,:,nZ) = zeros(nX,nY); % Because we can't resolve the last mode.

end


