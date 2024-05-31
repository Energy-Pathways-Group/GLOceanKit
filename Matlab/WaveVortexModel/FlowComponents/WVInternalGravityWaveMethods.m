classdef WVInternalGravityWaveMethods < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties (Abstract,GetAccess=public, SetAccess=public)
        Ap,Am,A0
    end
    properties (Abstract,GetAccess=public, SetAccess=protected)
        z
        UAp,VAp,WAp,NAp
        UAm,VAm,WAm,NAm
        ApmD,ApmN
        % Omega, iOmega
        Apm_TE_factor
    end
    methods (Abstract)
        addPrimaryFlowComponent(self,primaryFlowComponent)
    end
    properties (Dependent,GetAccess=public, SetAccess=protected)
        % returns the internal gravity wave flow component
        %
        % - Topic: Primary flow components
        % - Declaration: waveComponent
        % - Returns flowComponent: subclass of WVPrimaryFlowComponent
        % - nav_order: 2
        waveComponent
    end
    methods (Abstract)
        ratio = maxFw(self,kMode,lMode,j)
    end

    methods (Access=protected)
        function initializeInternalGravityWaveComponent(self)
            % After the WVStratifiedFlow and WVTransform constructors have
            % finishes, this should be called to finish initialization of
            % this flow component.
            arguments
                self WVTransform
            end
            flowComponent = WVInternalGravityWaveComponent(self);
            self.addPrimaryFlowComponent(flowComponent);

            function initVariable(varName,value)
                if isempty(self.(varName)) || isscalar(self.(varName))
                    self.(varName) = value;
                else
                    self.(varName) = self.(varName) + value;
                end

            end
            [ApmD_,ApmN_] = flowComponent.internalGravityWaveSpectralTransformCoefficients;
            [UAp_,VAp_,WAp_,NAp_] = flowComponent.internalGravityWaveSpatialTransformCoefficients;
            initVariable("ApmD",ApmD_);
            initVariable("ApmN",ApmN_);
            initVariable("UAp",UAp_);
            initVariable("VAp",VAp_);
            initVariable("WAp",WAp_);
            initVariable("NAp",NAp_);

            initVariable("UAm",conj(UAp_));
            initVariable("VAm",conj(VAp_));
            initVariable("WAm",WAp_);
            initVariable("NAm",-NAp_);

            self.iOmega = sqrt(-1)*self.Omega;

            initVariable("Apm_TE_factor",flowComponent.totalEnergyFactorForCoefficientMatrix(WVCoefficientMatrix.Ap));

            self.addVariableAnnotations(WVInternalGravityWaveMethods.variableAnnotationsForInternalGravityWaveComponent);
        end
    end

    methods
        function flowComponent = get.waveComponent(self)
            flowComponent = self.flowComponent('wave');
        end

        function energy = waveEnergy(self)
            % total energy of the geostrophic flow
            %
            % - Topic: Energetics
            % - Declaration: geostrophicEnergy
            % - nav_order: 2
            energy = self.totalEnergyOfFlowComponent(self.flowComponent('wave'));
        end

        function [omega,k,l] = initWithWaveModes(self, options)
            % initialize with the given wave modes
            %
            % $$
            % sin(k*x+l*y)*F_j*sin(omega*t + phi)
            % $$
            %
            % Clears variables Ap,Am,A0 and then sets the given wave modes.
            % - Topic: Initial conditions — Waves
            % - Declaration: [omega,k,l] = initWithWaveModes(kMode, lMode, j, phi, u, sign)
            % - Parameter kMode: integer index, (k0 > -Nx/2 && k0 < Nx/2)
            % - Parameter lMode: integer index, (l0 > -Ny/2 && l0 < Ny/2)
            % - Parameter j: integer index, (j0 >= 1 && j0 <= nModes), unless k=l=0 in which case j=0 is okay (inertial oscillations)
            % - Parameter phi: phase in radians, (0 <= phi <= 2*pi)
            % - Parameter u: fluid velocity (m/s)
            % - Parameter sign: sign of the frequency, +1 or -1
            % - Returns omega: frequencies of the waves (radians/s)
            % - Returns k: wavenumber k of the waves (radians/m)
            % - Returns l: wavenumber l of the waves (radians/m)
            arguments
                self WVTransform {mustBeNonempty}
                options.kMode (:,1) double
                options.lMode (:,1) double
                options.j (:,1) double
                options.phi (:,1) double = 0
                options.u (:,1) double
                options.sign (:,1) double
            end
            self.Ap = zeros(self.spectralMatrixSize);
            self.Am = zeros(self.spectralMatrixSize);
            self.A0 = zeros(self.spectralMatrixSize);

            [omega,k,l] = self.setWaveModes(kMode=options.kMode,lMode=options.lMode,j=options.j,phi=options.phi,u=options.u,sign=options.sign);
        end

        function [omega,k,l] = setWaveModes(self, options)
            % set amplitudes of the given wave modes
            %
            % Overwrite any existing wave modes with the given new values
            % - Topic: Initial conditions — Waves
            % - Declaration: [omega,k,l] = setWaveModes(kMode, lMode, j, phi, u, sign)
            % - Parameter kMode: integer index, (k0 > -Nx/2 && k0 < Nx/2)
            % - Parameter lMode: integer index, (l0 > -Ny/2 && l0 < Ny/2)
            % - Parameter j: integer index, (j0 >= 1 && j0 <= nModes), unless k=l=0 in which case j=0 is okay (inertial oscillations)
            % - Parameter phi: phase in radians, (0 <= phi <= 2*pi)
            % - Parameter u: fluid velocity (m/s)
            % - Parameter sign: sign of the frequency, +1 or -1
            % - Returns omega: frequencies of the waves (radians/s)
            % - Returns k: wavenumber k of the waves (radians/m)
            % - Returns l: wavenumber l of the waves (radians/m)
            arguments
                self WVTransform {mustBeNonempty}
                options.kMode (:,1) double
                options.lMode (:,1) double
                options.j (:,1) double
                options.phi (:,1) double = 0
                options.u (:,1) double
                options.sign (:,1) double
            end

            [kMode,lMode,j,u,phi,omegasign] = self.waveComponent.normalizeWaveModeProperties(options.kMode, options.lMode, options.j, options.u, options.phi, options.sign);
            indices = self.indexFromModeNumber(kMode,lMode,j);
            phi = phi - omegasign.*self.Omega(indices)*(self.t-self.t0);
            A = u.*exp(sqrt(-1)*phi)/(2*self.maxFw(kMode,lMode,j));

            Ap_ = self.Ap;
            Am_ = self.Am;
            Ap_(indices(omegasign>0)) = A(omegasign>0);
            Am_(indices(omegasign<0)) = A(omegasign<0);
            self.throwErrorIfDensityViolation(A0=self.A0,Ap=Ap_,Am=Am_,additionalErrorInfo=sprintf('The modes you are setting will cause the fluid state to violate this condition.\n'));
            self.Ap(indices(omegasign>0)) = A(omegasign>0);
            self.Am(indices(omegasign<0)) = A(omegasign<0);

            k = self.K(indices);
            l = self.L(indices);
            omega = self.Omega(indices);
        end

        function [omega,k,l] = addWaveModes(self, options)
            % add amplitudes of the given wave modes
            %
            % Add new amplitudes to any existing amplitudes
            % - Topic: Initial conditions — Waves
            % - Declaration: [omega,k,l] = addWaveModes(kMode, lMode, j, phi, u, sign)
            % - Parameter kMode: integer index, (k0 > -Nx/2 && k0 < Nx/2)
            % - Parameter lMode: integer index, (l0 > -Ny/2 && l0 < Ny/2)
            % - Parameter j: integer index, (j0 >= 1 && j0 <= nModes), unless k=l=0 in which case j=0 is okay (inertial oscillations)
            % - Parameter phi: phase in radians, (0 <= phi <= 2*pi)
            % - Parameter u: fluid velocity (m/s)
            % - Parameter sign: sign of the frequency, +1 or -1
            % - Returns omega: frequencies of the waves (radians/s)
            % - Returns k: wavenumber k of the waves (radians/m)
            % - Returns l: wavenumber l of the waves (radians/m)
            arguments
                self WVTransform {mustBeNonempty}
                options.kMode (:,1) double
                options.lMode (:,1) double
                options.j (:,1) double
                options.phi (:,1) double = 0
                options.u (:,1) double
                options.sign (:,1) double
            end

            [kMode,lMode,j,u,phi,omegasign] = self.waveComponent.normalizeWaveModeProperties(options.kMode, options.lMode, options.j, options.u, options.phi, options.sign);
            indices = self.indexFromModeNumber(kMode,lMode,j);
            phi = phi - omegasign.*self.Omega(indices)*(self.t-self.t0);
            A = u.*exp(sqrt(-1)*phi)/(2*self.maxFw(kMode,lMode,j));

            Ap_ = self.Ap;
            Am_ = self.Am;
            Ap_(indices(omegasign>0)) = Ap_(indices(omegasign>0)) + A(omegasign>0);
            Am_(indices(omegasign<0)) = Am_(indices(omegasign<0)) + A(omegasign<0);
            self.throwErrorIfDensityViolation(A0=self.A0,Ap=Ap_,Am=Am_,additionalErrorInfo=sprintf('The modes you are setting will cause the fluid state to violate this condition.\n'));
            self.Ap = Ap_;
            self.Am = Am_;

            k = self.K(indices);
            l = self.L(indices);
            omega = self.Omega(indices);
        end

        function removeAllWaves(self)
            % removes all wave from the model, including inertial oscillations
            %
            % Simply sets Ap and Am to zero.
            % - Topic: Initial conditions — Waves
            self.Ap(logical(self.waveComponent.maskAp)) = 0;
            self.Am(logical(self.waveComponent.maskAm)) = 0;
        end

        function [omega, alpha, k, l, mode, phi, A, norm] = waveModesFromWaveCoefficients(self)
            % Returns normalized amplitudes and phases of all waves
            %
            % This returns the properties of the waves being used in the
            % gridded simulation, as their properly normalized individual
            % wave components. Very useful for debugging.
            %
            % Note that A_plus and A_minus each have half the inertial
            % energy. This can be misleading, but the phasing is chosen to
            % make it work. Never-the-less, we double/zero that component.
            %
            % - Topic: Initial conditions — Waves
            % - Declaration: [omega, alpha, k, l, mode, phi, A, norm] = waveModesFromWaveCoefficients()
            A_p = self.Ap;
            A_p(1,1,:) = 2*A_p(1,1,:);
            A_m = self.Am;
            A_m(1,1,:) = 0*A_m(1,1,:);

            [K,L,J] = self.kljGrid;
            Omega = self.Omega;

            [A_plus,phi_plus,linearIndex] = WVTransform.extractNonzeroWaveProperties(A_p);
            omega_plus = Omega(linearIndex);
            mode_plus = J(linearIndex);
            alpha_plus = atan2(L(linearIndex),K(linearIndex));
            k_plus = K(linearIndex);
            l_plus = L(linearIndex);

            [A_minus,phi_minus,linearIndex] = WVTransform.extractNonzeroWaveProperties(A_m);
            omega_minus = -Omega(linearIndex);
            mode_minus = J(linearIndex);
            alpha_minus = atan2(L(linearIndex),K(linearIndex));
            k_minus = K(linearIndex);
            l_minus = L(linearIndex);

            k = [k_plus; k_minus];
            l = [l_plus; l_minus];
            omega = [omega_plus; omega_minus];
            mode = [mode_plus; mode_minus];
            alpha = [alpha_plus; alpha_minus];
            phi = [phi_plus; phi_minus];
            A = [A_plus; A_minus];
            norm = Normalization.kConstant;
        end
        function initWithGMSpectrum(self, GMAmplitude, varargin)
            % initialize with a Garrett-Munk spectrum
            %
            % This only initializes the wave components, A0 is left untouched.
            %
            % - Topic: Initial conditions — Waves
            % - Declaration: initWithGMSpectrum( GMAmplitude, varargin)
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
            B_norm = 1/atan(sqrt(self.Nmax*self.Nmax/(self.f*self.f)-1));

            % This function tells you how much energy you need between two
            % frequencies for a given vertical mode.
            GM2D_int = @(omega0,omega1,j) E*H_norm*B_norm*((j+j_star).^(-5/2))*(atan(self.f/sqrt(omega0*omega0-self.f*self.f)) - atan(self.f/sqrt(omega1*omega1-self.f*self.f)));

            % Do a quick check to see how much energy is lost due to
            % limited vertical resolution.
            maxMode = self.Nj;
            for iArg = 1:2:length(varargin)
                if strcmp(varargin{iArg}, 'maxMode')
                    maxMode = varargin{iArg+1};
                end
            end

            totalEnergy = 0;
            for mode=1:maxMode
                totalEnergy = totalEnergy + GM2D_int(self.f,self.Nmax,mode);
            end
            fprintf('You will miss %.2f%% of the energy due to limited vertical modes.\n',100-100*totalEnergy/E);

            [GM3Dint,GM3Dext] = self.initWithSpectralFunction(GM2D_int,varargin{:});

            fprintf('After distributing energy across frequency and mode, you still have %.2f%% of reference GM energy.\n',100*(sum(sum(sum(GM3Dint))) + sum(GM3Dext))/E);
            fprintf('Due to restricted domain size, the j=1,k=l=0 mode contains %.2f%% the total energy.\n',100*GM3Dint(1,1,1)/(sum(sum(sum(GM3Dint))) + sum(GM3Dext)) );

            GM_sum_int = sum(sum(sum(GM3Dint)))/E;
            GM_sum_ext = sum(GM3Dext)/E;
            C = self.Apm_TE_factor.*(self.Ap.*conj(self.Ap) + self.Am.*conj(self.Am));
            GM_random_sum_int = sum( C(:) )/E;
            GM_random_sum_ext = sum((self.offgridModes.U_cos_ext.*self.offgridModes.U_cos_ext + self.offgridModes.V_cos_ext.*self.offgridModes.V_cos_ext).*self.offgridModes.h_ext/2)/E;
            fprintf('The (gridded, external) wave field sums to (%.2f%%, %.2f%%) GM given the scales, and the randomized field sums to (%.2f%%, %.2f%%) GM\n', 100*GM_sum_int, 100*GM_sum_ext, 100*GM_random_sum_int,100*GM_random_sum_ext);
        end

        function [GM3Dint,GM3Dext] = initWithSpectralFunction(self, GM2D_int, varargin)
            % initialize the wave spectrum with a given function
            %
            % - Topic: Initial conditions — Waves
            % - Declaration: [GM3Dint,GM3Dext] = initWithSpectralFunction(GM2D_int, varargin)
            %
            % The GM2D_int function is used to assign variance to a given
            % wave mode. It has three arguments, omega0, omega1, and j and
            % should return the amount of variance you want assigned to a
            % wave mode between omega0 and omega1 at vertical mode j.
            %
            % The returned values GM3Dint are the results of distributing
            % this variance. size(GM3Dint) = size(Kh), so you can see
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
            % integrated over for assigned energy. By default it is self.Nmax-self.f
            if nargin(GM2D_int) ~= 3
                error('The spectral function must take three inputs: omega0, omega1, and j.\n');
            end

            if mod(length(varargin),2) ~= 0
                error('Arguments must be given as name/value pairs.');
            end

            [K,L,~] = self.kljGrid;
            Kh = sqrt(K.*K + L.*L);
            Omega = self.Omega;

            % Set defaults
            shouldRandomizeAmplitude = 0;
            maxDeltaOmega = self.Nmax-self.f;
            initializeModes = 0;
            energyWarningThreshold = 0.5;
            excludeNyquist = 1;
            minK = 0;
            maxK = max(max(max(abs(K))));
            minMode = 1;
            maxMode = self.Nj-1; % j=0 mode shouldn't be used/counted

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
                nyquistIndicesForK = sub2ind(size(Omega),repmat((ceil(self.Nx/2)+1)*ones(1,self.Ny),[1 self.Nj]),repmat(1:self.Ny,[1 self.Nj]),reshape(ones(1,self.Ny)'*(1:self.Nj),1,[]));
                nyquistIndicesForL = sub2ind(size(Omega),repmat(1:self.Nx,[1 self.Nj]),repmat((ceil(self.Ny/2)+1)*ones(1,self.Nx),[1 self.Nj]),reshape(ones(1,self.Nx)'*(1:self.Nj),1,[]));
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
            internalOmegaLinearIndices = reshape(1:numel(Omega),size(Omega));
            externalOmegaLinearIndices = 1:length(self.offgridModes.omega_ext);
            GM3Dint = zeros(size(Kh));
            GM3Dext = zeros(size(self.offgridModes.k_ext));
            indicesToSkip = reshape(internalOmegaLinearIndices(abs(Kh) < minK | abs(Kh) > maxK),[],1);
            for iMode = minMode:maxMode
                intOmegasLinearIndicesForIMode = reshape(internalOmegaLinearIndices(:,:,iMode+1),[],1); % Find the indices of everything with this iMode (iMode+1 b/c of j=0 mode)
                if excludeNyquist == 1
                    intOmegasLinearIndicesForIMode = setdiff( intOmegasLinearIndicesForIMode, nyquistIndices); % Now remove indices associated with the Nyquist
                end
                intOmegasLinearIndicesForIMode = setdiff( intOmegasLinearIndicesForIMode, indicesToSkip); % And remove indices outside the user specified wavenumber threshold.

                % Flatten the internal omegas (and their index)
                intOmegas = abs(Omega(intOmegasLinearIndicesForIMode));

                % Now do the same for the external modes
                indices = find(self.offgridModes.j_ext == iMode);
                extOmegas = reshape(abs(self.offgridModes.omega_ext(indices)),[],1);
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
                totalEnergyInThisMode = GM2D_int(self.f,self.Nmax,iMode);

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
                    leftDeltaOmega = omega0 - self.f;
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
            A = sqrt((GM3Dint./self.h)/2);

            % Each standard coefficient a(i,j,k) has equal conjugate, which
            % is already accounted for as having the same frequency. So
            % sum(a^2/2) really is the total energy.
            if shouldRandomizeAmplitude == 1
                A_plus = A.*self.generateHermitianRandomMatrix( );
                A_minus = A.*self.generateHermitianRandomMatrix(  );

                self.offgridModes.U_ext = sqrt(2*GM3Dext./self.offgridModes.h_ext).*randn( size(self.offgridModes.h_ext) );
                self.offgridModes.PrecomputeExternalWaveCoefficients();
            else
                % Randomize phases, but keep unit length
                A_plus = self.generateHermitianRandomMatrix(  shouldExcludeNyquist=excludeNyquist, allowMeanPhase=1 );
                A_minus = self.generateHermitianRandomMatrix( shouldExcludeNyquist=excludeNyquist, allowMeanPhase=1 );

                goodIndices = abs(A_plus) > 0;
                A_plus(goodIndices) = A_plus(goodIndices)./abs(A_plus(goodIndices));
                A_plus = A.*A_plus;
                goodIndices = abs(A_minus) > 0;
                A_minus(goodIndices) = A_minus(goodIndices)./abs(A_minus(goodIndices));
                A_minus = A.*A_minus;

                % Check this factor of 2!!! Is the correct? squared
                % velocity to energy, I think.
                self.offgridModes.U_ext = sqrt(2*GM3Dext./self.offgridModes.h_ext);
                self.offgridModes.PrecomputeExternalWaveCoefficients();
            end

            A_minus(1,1,:) = conj(A_plus(1,1,:)); % Inertial motions go only one direction!
            self.Ap = A_plus;
            self.Am = A_minus;
        end

        function initWithHorizontalWaveNumberSpectrum(self,GMAmplitude,options)
            % initialize with a Alternative Interal Wave Spectrum in
            % function of horizontal wave number and mode
            %
            % This only initializes the wave components, A0 is left untouched.
            %
            % - Topic: Initial conditions — Waves
            % - Declaration: initWithHorizontalWaveNumberSpectrum(GMAmplitude,options)
            % - Parameter GMAmplitude:
            % - Parameter j_star: (optional)
            % - Parameter slope: (optional)


            arguments
                self WVTransform {mustBeNonempty}
                GMAmplitude (1,1) double
                options.j_star (1,1) double = 3
                options.slope (1,1) double = 1
                options.shouldRandomizeAmplitude = 1
            end


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Create a reasonable total wavenumber (Radial) axis

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Kh= self.Kh;
            kRadial = self.radialWavenumberAxis();
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Distribution of Energy

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            j_star=options.j_star;
            slope=options.slope;

            % GM Parameters
            L_gm = 1.3e3; % thermocline exponential scale, meters
            invT_gm = 5.2e-3; % reference buoyancy frequency, radians/seconds
            E_gm = 6.3e-5; % non-dimensional energy parameter
            E_T = L_gm*L_gm*L_gm*invT_gm*invT_gm*E_gm*GMAmplitude;
            %  E = E*(self.Lz/L_gm); % This correction fixes the amplitude so that the HKE variance at a given depth matches (instead of depth integrated energy)

            % Compute the proper vertical function normalization
            M = (j_star^2 +(2:1024).^2).^((-5/4));
            M_norm = sum(M);

            %Create the energy matrix 3D
            TotalEnergy = zeros(size(Kh));
            Energy_slice=zeros(size(Kh(:,:,1)));

            %%%% Redistributing the energy %%%
            for j=(2:length(self.j))
                h=self.h(1,1,j);
                Kh_2D = Kh(:,:,j);
                LR= sqrt(self.g*h)/self.f;

                fun = @(k) (1./(k.^2*LR^2 + 1).^(1*slope))*LR;
                B_norm = integral(fun,kRadial(1),kRadial(end));

                for i=(1:length(kRadial)-1)

                    % Integrate the energy btw 2 Kh
                    E = E_T*(integral(fun,kRadial(i),kRadial(i+1))/B_norm)*(((j^2 + j_star^2).^((-5/4)))/M_norm);

                    % find all the kl point btw the two values of Kh
                    ind = find(Kh_2D>=kRadial(i) & Kh_2D<kRadial(i+1));

                    % Distribuite equally the energy btw all the points that are
                    % btw the circles (values of Kh)
                    n = length(ind);
                    if n > 0
                        Energy_slice(ind) = E/n;
                    end
                end
                TotalEnergy(:,:,j)= Energy_slice;

                %clear the slice but maintain the size
                Energy_slice=zeros(size(Kh(:,:,1)));
            end


            disp([' Initial Total energy:',num2str(E_T), ', Total energy after distribution:',num2str(sum(TotalEnergy(:)))])
            fprintf('After distributing energy across frequency and mode, you still have %.2f%% of reference GM energy.\n',100*sum(TotalEnergy(:))/E_T);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % After comput the amplitude I insert that in the model
            % to get the variables in real space

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %The amplitude is:

            shouldRandomizeAmplitude=options.shouldRandomizeAmplitude;

            A = sqrt((TotalEnergy./self.h)/2);

            if shouldRandomizeAmplitude == 1
                A_plus = A.*self.generateHermitianRandomMatrix( );
                A_minus = A.*self.generateHermitianRandomMatrix(  );

                %self.offgridModes.U_ext = sqrt(2*GM3Dext./self.offgridModes.h_ext).*randn( size(self.offgridModes.h_ext) );
                %self.offgridModes.PrecomputeExternalWaveCoefficients();
            else
                % Randomize phases, but keep unit length
                A_plus = self.generateHermitianRandomMatrix(shouldExcludeNyquist=1, allowMeanPhase=1 );
                A_minus = self.generateHermitianRandomMatrix(shouldExcludeNyquist=1, allowMeanPhase=1 );

                goodIndices = abs(A_plus) > 0;
                A_plus(goodIndices) = A_plus(goodIndices)./abs(A_plus(goodIndices));
                A_plus = A.*A_plus;

                goodIndices = abs(A_minus) > 0;
                A_minus(goodIndices) = A_minus(goodIndices)./abs(A_minus(goodIndices));
                A_minus = A.*A_minus;

                % Check this factor of 2!!! Is the correct? squared
                % velocity to energy, I think.
                %self.offgridModes.U_ext = sqrt(2*GM3Dext./self.offgridModes.h_ext);
                %self.offgridModes.PrecomputeExternalWaveCoefficients();
            end

            self.Ap=A_plus;
            self.Am=A_minus;

        end
    end
    methods (Static, Hidden=true)
        function variableAnnotations = variableAnnotationsForInternalGravityWaveComponent()
            % return array of WVVariableAnnotation instances initialized by default
            %
            % This function creates annotations for the built-in variables supported by
            % the WVTransform.
            %
            % - Topic: Internal
            % - Declaration: operations = defaultVariableAnnotations()
            % - Returns operations: array of WVVariableAnnotation instances

            variableAnnotations = WVVariableAnnotation.empty(0,0);

            annotation = WVVariableAnnotation('waveEnergy',{},'m3/s2', 'total energy, waves');
            annotation.isVariableWithLinearTimeStep = 0;
            annotation.isVariableWithNonlinearTimeStep = 1;
            variableAnnotations(end+1) = annotation;
        end

        function [A,phi,linearIndex] = extractNonzeroWaveProperties(Matrix)
            % Takes a Hermitian matrix and returns the amplitude and phase of nonzero components
            %
            % This function makes assumptions about the structure of the matrix.
            % - Topic: Utility function
            % - Declaration: [A,phi,linearIndex] = ExtractNonzeroWaveProperties(Matrix)
            % - Parameter Matrix: Hermitian conjugate matrix
            % - Returns A: amplitude
            % - Returns phi: phase
            % - Returns linearIndex: linear index of matrix component
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

    end
end