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
        
        K2, Kh, h, Omega, f0, C, rho0
        u_plus, u_minus, v_plus, v_minus, w_plus, w_minus, zeta_plus, zeta_minus
        
        Amp_plus, Amp_minus % Used for diagnostics only
        
        U_cos_int, U_sin_int, V_cos_int, V_sin_int, W_sin_int, Zeta_cos_int, Rho_cos_int, k_int, l_int, j_int, h_int, omega_int, phi_int

        U_cos_ext, U_sin_ext, V_cos_ext, V_sin_ext, W_sin_ext, Zeta_cos_ext, Rho_cos_ext, k_ext, l_ext, j_ext, k_z_ext, h_ext, omega_ext, phi_ext, F_ext, G_ext, norm_ext
        
        Xc, Yc, Zc % These may contain 'circular' versions of the grid
                
        version = 1.5
        performSanityChecks = 0
        advectionSanityCheck = 0
    end
    
    methods (Abstract, Access = protected)
        [F,G,h] = ModesAtWavenumber(self, k, norm ) % Return the normal modes and eigenvalue at a given wavenumber.
        [F,G,h] = ModesAtFrequency(self, omega, norm ) % Return the normal modes and eigenvalue at a given frequency.
        F = InternalUVModesAtDepth(self, z) % Returns normal modes at requested depth, size(F) = [length(z) nIntModes]
        G = InternalWModesAtDepth(self, z) % Returns normal modes at requested depth, size(G) = [length(z) nIntModes]
        F = ExternalUVModesAtDepth(self, z) % Returns normal modes at requested depth, size(F) = [length(z) nExtModes]
        G = ExternalWModesAtDepth(self, z) % Returns normal modes at requested depth, size(G) = [length(z) nExtModes]
        rho = RhoBarAtDepth(self,z)
        N2 = N2AtDepth(self,z)
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
            
            % Spectral domain, in radians
            dk = 1/self.Lx;          % fourier frequency
            self.k = 2*pi*([0:ceil(self.Nx/2)-1 -floor(self.Nx/2):-1]*dk)';
            
            dl = 1/self.Ly;          % fourier frequency
            self.l = 2*pi*([0:ceil(self.Ny/2)-1 -floor(self.Ny/2):-1]*dl)';
            
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
        
        function [omega, alpha, mode, phi, A] = WaveCoefficientsFromGriddedWaves(self)
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
            
            [A_minus,phi_minus,linearIndex] = ExtractNonzeroWaveProperties(A_m);
            omega_minus = -self.Omega(linearIndex);
            mode_minus = self.J(linearIndex);
            alpha_minus = atan2(self.L(linearIndex),self.K(linearIndex));
            
            omega = [omega_plus; omega_minus];
            mode = [mode_plus; mode_minus];
            alpha = [alpha_plus; alpha_minus];
            phi = [phi_plus; phi_minus];
            A = [A_plus; A_minus];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % These functions add free waves to the model, not constrained to
        % the FFT grid.
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function omega = SetExternalWavesWithWavenumbers(self, k, l, j, phi, A, type)
            % Add free waves to the model that are not constrained to the
            % gridded solution. The amplitude, A, can be given as an energy
            % density or a maximum U velocity. The type must then be either
            % 'energyDensity' or 'maxU'.            
            self.k_ext = reshape(k,1,[]);
            self.l_ext = reshape(l,1,[]);
            self.j_ext = reshape(j,1,[]);
            self.k_z_ext = self.j_ext*pi/self.Lz;
            self.phi_ext = reshape(phi,1,[]);
                     
            if strcmp(type, 'energyDensity')
                self.norm_ext = 'const_G_norm';
            elseif strcmp(type, 'maxU')
                self.norm_ext = 'max_u';
            end
            
            
            K2h = self.k_ext.*self.k_ext + self.l_ext.*self.l_ext;
            
            self.h_ext = zeros(size(self.j_ext));
            self.F_ext = zeros(length(self.z),length(self.j_ext));
            self.G_ext = zeros(length(self.z),length(self.j_ext));
            for iWave=1:length(j)
                [FExt,GExt,hExt] = self.ModesAtWavenumber(sqrt(K2h(iWave)), self.norm_ext);
                self.F_ext(:,iWave) = FExt(:,j(iWave));
                self.G_ext(:,iWave) = GExt(:,j(iWave));
                self.h_ext(iWave) = hExt(j(iWave));
            end
            
            U = reshape(A,1,[]);
            if strcmp(type, 'energyDensity')
                U = U./sqrt(self.h_ext);
            end
            
            g = 9.81;
            self.omega_ext = sqrt( g*self.h_ext .* K2h + self.f0*self.f0 );
            omega = self.omega_ext;
            
            self.PrecomputeExternalWaveCoefficients(U);           
        end
        
        function k = SetExternalWavesWithFrequencies(self, omega, alpha, j, phi, A, type)
            % Add free waves to the model that are not constrained to the
            % gridded solution. The amplitude, A, can be given as an energy
            % density or a maximum U velocity. The type must then be either
            % 'energyDensity' or 'maxU'.                       
            self.j_ext = reshape(j,1,[]);
            self.k_z_ext = self.j_ext*pi/self.Lz;
            self.omega_ext = reshape(omega,1,[]);
            self.phi_ext = reshape(phi,1,[]);

            if strcmp(type, 'energyDensity')
                self.norm_ext = 'const_G_norm';
            elseif strcmp(type, 'maxU')
                self.norm_ext = 'max_u';
            end
                        
            self.h_ext = zeros(size(self.j_ext));
            self.F_ext = zeros(length(self.z),length(self.j_ext));
            self.G_ext = zeros(length(self.z),length(self.j_ext));
            for iWave=1:length(j)
                [FExt,GExt,hExt] = self.ModesAtFrequency(omega(iWave), self.norm_ext);
                self.F_ext(:,iWave) = FExt(:,j(iWave));
                self.G_ext(:,iWave) = GExt(:,j(iWave));
                self.h_ext(iWave) = hExt(j(iWave));
            end
            
            U = reshape(A,1,[]);
            if strcmp(type, 'energyDensity')
                U = U./sqrt(self.h_ext);
            end
            
            g = 9.81;
            k = sqrt((self.omega_ext.*self.omega_ext - self.f0*self.f0)./(g*self.h_ext));
            alpha0 = reshape(alpha,1,[]);
            self.k_ext = k .* cos(alpha0);
            self.l_ext = k .* sin(alpha0);
            
            self.PrecomputeExternalWaveCoefficients(U);
        end
        
        function PrecomputeExternalWaveCoefficients(self, U)
            alpha0 = atan2(self.l_ext,self.k_ext);
            Kh_ = sqrt( self.k_ext.*self.k_ext + self.l_ext.*self.l_ext);
            
            self.U_cos_ext = U .* cos(alpha0);
            self.U_sin_ext = U .* (self.f0 ./ self.omega_ext) .* sin(alpha0);
            self.V_cos_ext = U .* sin(alpha0);
            self.V_sin_ext = -U .* (self.f0 ./ self.omega_ext) .* cos(alpha0);
            self.W_sin_ext = U .* Kh_ .* self.h_ext;
            self.Zeta_cos_ext = - U .* Kh_ .* self.h_ext ./ self.omega_ext;
            self.Rho_cos_ext = - (self.rho0/9.81) .* U .* Kh_ .* self.h_ext ./ self.omega_ext;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Create a full Garrett-Munk spectrum (public)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function InitializeWithGMSpectrum(self, Amp, shouldRandomizeAmplitude)
            if ~exist('shouldRandomizeAmplitude', 'var')
                shouldRandomizeAmplitude = 0;
            end
            
            % GM Parameters
            j_star = 3;
            L_gm = 1.3e3; % thermocline exponential scale, meters
            invT_gm = 5.2e-3; % reference buoyancy frequency, radians/seconds
            E_gm = 6.3e-5; % non-dimensional energy parameter
            E = L_gm*L_gm*L_gm*invT_gm*invT_gm*E_gm*Amp;
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
            totalEnergy = 0;
            for mode=1:self.nModes
                totalEnergy = totalEnergy + GM2D_int(self.f0,self.N0,mode);
            end
            fprintf('You are missing %.2f%% of the energy due to limited vertical modes.\n',100-100*totalEnergy/E);
                        
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
            for iMode = 1:self.nModes
                % Flatten the internal omegas (and their index)
                intOmegas = reshape(abs(self.Omega(:,:,iMode)),[],1);
                intOmegasLinearIndicesForIMode = reshape(internalOmegaLinearIndices(:,:,iMode),[],1);
                
                % Now do the same for the external modes
                indices = find(self.j_ext == iMode);
                extOmegas = reshape(abs(self.omega_ext(indices)),[],1);
                extOmegasLinearIndicesForIMode = reshape(externalOmegaLinearIndices(indices),[],1);
                
                % Make a combined list, but note which list each omega came
                % from.
                allOmegas = cat(1,intOmegas,extOmegas);
                allIndices = cat(1,intOmegasLinearIndicesForIMode,extOmegasLinearIndicesForIMode);
                allSource = cat(1,zeros(size(intOmegas)), ones(size(extOmegas)));
                            
                % Sort the frequencies for this mode.
                [sortedOmegas, sortedOmegasIndices] = sort(allOmegas);
                sortedIndices = allIndices(sortedOmegasIndices);
                sortedSource = allSource(sortedOmegasIndices);
                
                % Then find where the omegas differ.
                omegaDiffIndices = find(diff(sortedOmegas) > 0);               
            
                lastIdx = 1;
                omega0 = sortedOmegas(lastIdx);
                leftDeltaOmega = 0;
                for idx=omegaDiffIndices'
                    currentIdx = idx+1;
                    nOmegas = currentIdx-lastIdx;
                                 
                    omega1 = sortedOmegas(idx + 1);
                    rightDeltaOmega = (omega1-omega0)/2;
                    energyPerFrequency = GM2D_int(omega0-leftDeltaOmega,omega0+rightDeltaOmega,iMode)/nOmegas;
                    
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
            end
            fprintf('After distributing energy across frequency and mode, you still have %.2f%% of reference GM energy.\n',100*(sum(sum(sum(GM3Dint))) + sum(GM3Dext))/E);
            fprintf('Due to restricted domain size, the j=1,k=l=0 mode contains %.2f%% the total energy.\n',100*GM3Dint(1,1,1)/(sum(sum(sum(GM3Dint))) + sum(GM3Dext)) );
            
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
                
                self.phi_ext = 2*pi*rand( size(self.k_ext) );
                U = sqrt(2*GM3Dext/self.hExternal).*randn( size(self.uExternal) );
                self.PrecomputeExternalWaveCoefficients(U);                
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
                self.phi_ext = 2*pi*rand( size(self.k_ext) );
                U = sqrt(2*GM3Dext./self.h_ext);
                self.PrecomputeExternalWaveCoefficients(U);   
            end
            
            
            A_minus(1,1,:) = conj(A_plus(1,1,:)); % Intertial motions go only one direction!
            
            GM_sum_int = sum(sum(sum(GM3Dint)))/E;
            GM_sum_ext = sum(GM3Dext)/E;
            GM_random_sum_int = sum(sum(sum(A_plus.*conj(A_plus) + A_minus.*conj(A_minus)  )))/E;
            GM_random_sum_ext = sum((self.U_cos_ext.*self.U_cos_ext + self.V_cos_ext.*self.V_cos_ext).*self.h_ext/2)/E;
            fprintf('The (gridded, external) wave field sums to (%.2f%%, %.2f%%) GM given the scales, and the randomized field sums to (%.2f%%, %.2f%%) GM\n', 100*GM_sum_int, 100*GM_sum_ext, 100*GM_random_sum_int,100*GM_random_sum_ext);
            
            self.GenerateWavePhases(A_plus,A_minus);
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
            
            self.SetExternalWavesWithFrequencies(omegaExt,alphaExt,jExt,zeros(size(omegaExt)),zeros(size(omegaExt)),'energyDensity');
            
            fprintf('Added %d external waves to fill out the GM spectrum.\n', length(omegaExt));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Return the dynamical fields on the grid at a given time
        % (Eulerian)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [u,v,w] = VelocityFieldAtTime(self, t)
            % Return the velocity field, which is the sum of the gridded
            % and external/free waves at time t. Note that if you do not
            % need w, don't request it and it won't be computed.
            
            % Get the gridded velocity field...
            if nargout == 3
                [u, v, w] = self.InternalVelocityFieldAtTime(t);
            else
                [u, v] = self.InternalVelocityFieldAtTime(t);
            end
            
            % ...add the external waves to the gridded solution.
            if ~isempty(self.k_ext)
                if nargout == 3
                    [u_ext, v_ext, w_ext] = self.ExternalVelocityFieldsAtTime(t);
                    w = w + w_ext;
                else
                    [u_ext, v_ext] = self.ExternalVelocityFieldsAtTime(t);
                end
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
            
            if ~isempty(self.k_ext)
                [w_ext, zeta_ext] = self.ExternalVerticalFieldsAtTime(t);
                w = w + w_ext;
                zeta = zeta + zeta_ext;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Return the dynamical fields at a given location and time
        % (Lagrangian)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [u] = VelocityAtTimePositionVector(self,t,p, varargin)
            % Return the velocity at time t and position p, where size(p) =
            % [3 n]. The rows of p represent [x,y,z].
            % Optional argument is passed to the interpolation method. I
            % recommend either linear or spline.
            psize = size(p);
            if psize(1) == 3
                [u,v,w] = self.VelocityAtTimePosition(t,p(1,:),p(2,:),p(3,:), varargin{:});
                u = cat(1,u,v,w);
            else
                [u,v,w] = self.VelocityAtTimePosition(t,p(:,1),p(:,2),p(:,3), varargin{:});
                u = cat(2,u,v,w);
            end
        end
        
        function [u] = DrifterVelocityAtTimePositionVector(self,t,p, varargin)
            % Return the velocity at time t and position p, where size(p) =
            % [3 n]. The rows of p represent [x,y,z].
            psize = size(p);
            if psize(1) == 3
                [u,v,~] = self.VelocityAtTimePosition(t,p(1,:),p(2,:),p(3,:), varargin{:});
                u = cat(1,u,v,zeros(1,psize(2)));
            else
                [u,v,~] = self.VelocityAtTimePosition(t,p(:,1),p(:,2),p(:,3), varargin{:});
                u = cat(2,u,v,zeros(psize(1),1));
            end
        end
               
        function [u,v,w] = VelocityAtTimePosition(self,t,x,y,z,varargin)
            if self.advectionSanityCheck == 0
               self.advectionSanityCheck = 1;
               if (self.z(end)-self.z(1)) ~= self.Lz
                   warning('Vertical domain does not span the full depth of the ocean. This will lead to NaNs when advected particles leave the resolved domain.')
               end
            end
            
            if nargin == 5
                method = 'spline';
            else
                if strcmp(varargin{1},'exact')
                    if nargout == 3
                        [u,v,w] = self.ExactVelocityAtTimePosition(t,x,y,z);
                    else
                        [u,v] = self.ExactVelocityAtTimePosition(t,x,y,z);
                    end
                    return
                else
                    method = varargin{1};
                end
            end
            
            % Return the velocity at time t, and positions x,y,z.
            [U,V,W] = self.InternalVelocityFieldAtTime(t);
            
            % (x,y) are periodic for the gridded solution
            x_tilde = mod(x,self.Lx);
            y_tilde = mod(y,self.Ly);
            
            u = interpn(self.X,self.Y,self.Z,U,x_tilde,y_tilde,z,method);
            v = interpn(self.X,self.Y,self.Z,V,x_tilde,y_tilde,z,method);
            if nargout == 3
                w = interpn(self.X,self.Y,self.Z,W,x_tilde,y_tilde,z,method);
                badParticles = isnan(u)|isnan(v)|isnan(w);
            else
                badParticles = isnan(u)|isnan(v);
            end
            
            % The above will fail (and return nan) when a particle is
            % between the last point and the first point.
            if any(badParticles)
                if isempty(self.Xc)
                    [self.Xc,self.Yc,self.Zc] = ndgrid([self.x;self.Lx],[self.y;self.Ly],self.z);
                end
                
                makeperiodic = @(A) cat(2,cat(1,A,A(1,:,:)), cat(1,A(:,1,:),A(1,1,:)));
                
                Uc = makeperiodic(U);
                Vc = makeperiodic(V);
                u(badParticles) = interpn(self.Xc,self.Yc,self.Zc,Uc,x_tilde(badParticles),y_tilde(badParticles),z(badParticles));
                v(badParticles) = interpn(self.Xc,self.Yc,self.Zc,Vc,x_tilde(badParticles),y_tilde(badParticles),z(badParticles));
                if nargout == 3
                    Wc = makeperiodic(W);
                    w(badParticles) = interpn(self.Xc,self.Yc,self.Zc,Wc,x_tilde(badParticles),y_tilde(badParticles),z(badParticles));
                end
            end
            
            if nargout == 3
                [u_ext,v_ext,w_ext] = self.ExternalVelocityAtTimePosition(t,x,y,z);
                u = u+u_ext;
                v = v+v_ext;
                w = w+w_ext;
            else
                [u_ext,v_ext] = self.ExternalVelocityAtTimePosition(t,x,y,z);
                u = u+u_ext;
                v = v+v_ext;
            end
        end
        
        function zeta = ZetaAtTimePosition(self,t,x,y,z)
            [zeta] = self.InternalZetaAtTimePositionExact(t,x,y,z);
            [~,zeta_ext] = self.ExternalVerticalFieldsAtTimePosition(t,x,y,z);
            
            zeta = zeta+zeta_ext;
        end
        
        function rho = DensityAtTimePosition(self,t,x,y,z)
            rho_bar = self.RhoBarAtDepth(z);
            rho_int = self.InternalDensityPertubationAtTimePositionExact(t,x,y,z);
            rho_ext = self.ExternalDensityPerturbationAtTimePosition(t,x,y,z);
            
            rho = rho_bar + rho_int + rho_ext;
        end
                
        function ShowDiagnostics(self)
            % Display various diagnostics about the simulation.
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
            
            self.GenerateWavePhasesForSpectralAdvection();
        end
    end
    
    methods (Access = protected)
        
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
        end
        
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
               
        function [u,v,w] = InternalVelocityFieldAtTime(self, t)
            % Returns the velocity field from the gridded solution/waves.
            phase_plus = exp(sqrt(-1)*self.Omega*t);
            phase_minus = exp(-sqrt(-1)*self.Omega*t);
            u_bar = self.u_plus.*phase_plus + self.u_minus.*phase_minus;
            v_bar = self.v_plus.*phase_plus + self.v_minus.*phase_minus;
            
            if self.performSanityChecks == 1
                CheckHermitian(u_bar);CheckHermitian(v_bar);
            end
            
            u = self.TransformToSpatialDomainWithF(u_bar);
            v = self.TransformToSpatialDomainWithF(v_bar);
            
            if nargout == 3
                w_bar = self.w_plus.*phase_plus + self.w_minus.*phase_minus;
                w = self.TransformToSpatialDomainWithG(w_bar);
            end
        end
                
        function [u,v,w] = ExternalVelocityFieldsAtTime(self, t)
            [u,v,w] = self.ExternalVelocityAtTimePosition(t,reshape(self.X,[],1),reshape(self.Y,[],1), reshape(self.Z,[],1));
            u = reshape(u,self.Nx,self.Ny,self.Nz);
            v = reshape(v,self.Nx,self.Ny,self.Nz);
            w = reshape(w,self.Nx,self.Ny,self.Nz);
        end
        
        function [w,zeta] = ExternalVerticalFieldsAtTime(self, t)       
            [w,zeta] = self.ExternalVerticalFieldsAtTimePosition(t,reshape(self.X,[],1),reshape(self.Y,[],1), reshape(self.Z,[],1));
            w = reshape(w,self.Nx,self.Ny,self.Nz);
            zeta = reshape(zeta,self.Nx,self.Ny,self.Nz);
        end
        
        function [u,v,w] = ExternalVelocityAtTimePosition(self,t,x,y,z)
            % Return the velocity field associated with the manually added
            % free waves at specified positions.
            if isempty(self.k_ext)
                u = zeros(size(x)); v=u; w=u; return;
            end
            theta = x * self.k_ext + y * self.l_ext + (self.omega_ext*t + self.phi_ext); % [N M] + [1 M]
            cos_theta = cos(theta); % [N M]
            sin_theta = sin(theta); % [N M]
            F = self.ExternalUVModesAtDepth(z); % [N M]
            
            u = sum( (self.U_cos_ext .* cos_theta + self.U_sin_ext .* sin_theta) .* F, 2);
            v = sum( (self.V_cos_ext .* cos_theta + self.V_sin_ext .* sin_theta) .* F, 2);
                        
            if nargout == 3
                G = self.ExternalWModesAtDepth(z); % [N M]
                w = sum( (self.W_sin_ext .* sin_theta) .* G, 2);
            end     
        end
   
        function [w,zeta] = ExternalVerticalFieldsAtTimePosition(self,t,x,y,z)
            if isempty(self.k_ext)
                w = zeros(size(x)); zeta=w; return;
            end
            theta = x * self.k_ext + y * self.l_ext + (self.omega_ext*t + self.phi_ext); % [N M] + [1 M]
            cos_theta = cos(theta); % [N M]
            sin_theta = sin(theta); % [N M]
            G = self.ExternalWModesAtDepth(z); % [N M]
            w = sum( (self.W_sin_ext .* sin_theta) .* G, 2);
            zeta = sum( (self.Zeta_cos_ext .* cos_theta) .* G, 2);
        end
        
        function rho = ExternalDensityPerturbationAtTimePosition(self,t,x,y,z)
            if isempty(self.k_ext)
                rho =  zeros(size(x)); return;
            end
            theta = x * self.k_ext + y * self.l_ext + (self.omega_ext*t + self.phi_ext); % [N M] + [1 M]
            cos_theta = cos(theta); % [N M]
            G = self.ExternalWModesAtDepth(z); % [N M]
            N2_ = self.N2AtDepth(z); % [N 1]
            
            rho = N2_ .* sum( (self.Rho_cos_ext .* cos_theta) .* G, 2);
        end
        
        % size(x) = [N 1]
        % size(phi) = [1 M]
        function [u,v,w] = InternalVelocityAtTimePositionExact(self,t,x,y,z)
            % Return the velocity field associated with the gridded
            % velocity field, but using spectral interpolation, rather than
            % the FFT grid.  
            if isempty(self.k_int)
                u = zeros(size(x)); v=u; w=u; return;
            end
 
            theta = x * self.k_int + y * self.l_int + (self.omega_int*t + self.phi_int); % [N M] + [1 M]
            cos_theta = cos(theta); % [N M]
            sin_theta = sin(theta); % [N M]
            F = self.InternalUVModesAtDepth(z); % [N M]
            
            u = sum( (self.U_cos_int .* cos_theta + self.U_sin_int .* sin_theta) .* F, 2);
            v = sum( (self.V_cos_int .* cos_theta + self.V_sin_int .* sin_theta) .* F, 2);
                        
            if nargout == 3
                G = self.InternalWModesAtDepth(z); % [N M]
                w = sum( (self.W_sin_int .* sin_theta) .* G, 2);
            end     
        end
        
        function zeta = InternalZetaAtTimePositionExact(self,t,x,y,z)
            if isempty(self.k_int)
                zeta =  zeros(size(x)); return;
            end
            
            theta = x * self.k_int + y * self.l_int + (self.omega_int*t + self.phi_int); % [N M] + [1 M]
            cos_theta = cos(theta); % [N M]
            G = self.InternalWModesAtDepth(z); % [N M]
            
            zeta = sum( (self.Zeta_cos_int .* cos_theta) .* G, 2);
        end
        
        function rho = InternalDensityPertubationAtTimePositionExact(self,t,x,y,z)
            if isempty(self.k_int)
                rho =  zeros(size(x)); return;
            end
            
            theta = x * self.k_int + y * self.l_int + (self.omega_int*t + self.phi_int); % [N M] + [1 M]
            cos_theta = cos(theta); % [N M]
            G = self.InternalWModesAtDepth(z); % [N M]
            N2_ = self.N2AtDepth(z); % [N 1]
            
            rho = N2_ .* sum( (self.Rho_cos_int .* cos_theta) .* G, 2);
        end
                
        function [u,v,w] = ExactVelocityAtTimePosition(self,t,x,y,z)
            % Uses coeffs to add up sines and cosines for an exact
            % solution.
            if nargout == 3
                [u,v,w] = self.InternalVelocityAtTimePositionExact(t,x,y,z);
                [u_ext,v_ext,w_ext] = self.ExternalVelocityAtTimePosition(t,x,y,z);
                
                u = u+u_ext;
                v = v+v_ext;
                w = w+w_ext;
            else
                [u,v] = self.InternalVelocityAtTimePositionExact(t,x,y,z);
                [u_ext,v_ext] = self.ExternalVelocityAtTimePosition(t,x,y,z);
                
                u = u+u_ext;
                v = v+v_ext;
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
A = MakeHermitian(randn(size) + sqrt(-1)*randn(size) )/sqrt(2);
% A(1,1,:) = 2*real(A(1,1,:)); % Double the zero frequency
A(1,1,:) = 2*A(1,1,:); % Double the zero frequency
A(nX/2+1,1,:) = -2*real(A(nX/2+1,1,:)); % Double the Nyquist frequency
A(1,nY/2+1,:) = -2*real(A(1,nY/2+1,:)); % Double the Nyquist frequency
A(nX/2+1,nY/2+1,:) = -2*real(A(nX/2+1,nY/2+1,:)); % Double the Nyquist frequency
A(:,:,nZ) = zeros(nX,nY); % Because we can't resolve the last mode.

end


