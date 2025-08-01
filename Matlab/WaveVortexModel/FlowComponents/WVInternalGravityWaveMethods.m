classdef WVInternalGravityWaveMethods < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties %(GetAccess=public, SetAccess=protected)
        UAp,VAp,WAp,NAp
        UAm,VAm,WAm,NAm
        ApmD,ApmN
        % Omega, iOmega

        iOmega
    end

    properties (Dependent,GetAccess=public, SetAccess=protected)
        % returns the internal gravity wave flow component
        %
        % - Topic: Primary flow components
        % - Declaration: waveComponent
        % - Returns flowComponent: subclass of WVPrimaryFlowComponent
        % - nav_order: 2
        waveComponent

        Omega
    end
    properties (Abstract)
        % z
        Ap,Am
        h_pm
    end
    methods (Abstract)        
        removeAll(self)
        ratio = maxFw(self,kMode,lMode,j)
        addPrimaryFlowComponent(self,primaryFlowComponent)
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

            % self.addOperation(self.operationForDynamicalVariable('u','v','w','eta','p',flowComponent=self.waveComponent));
        end
    end

    methods
        function flowComponent = get.waveComponent(self)
            flowComponent = self.flowComponentWithName('wave');
        end

        function [Apt,Amt] = waveCoefficientsAtTimeT(self)
            phase = exp(self.iOmega*(self.t-self.t0));
            Apt = self.Ap .* phase;
            Amt = self.Am .* conj(phase);
        end

        function Omega = get.Omega(self)
            Omega = sqrt(self.g*self.h_pm.*(self.K .* self.K + self.L .* self.L) + self.f*self.f);
        end

        function energy = waveEnergy(self)
            % total energy of the geostrophic flow
            %
            % - Topic: Energetics
            % - Declaration: geostrophicEnergy
            % - nav_order: 2
            energy = self.totalEnergyOfFlowComponent(self.flowComponentWithName('wave'));
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
            self.removeAll();

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
        
        function initWavesWithFrequencySpectrum(self,options)
            % initialize with waves of a specified frequency spectrum
            %
            % This allows you to initialize the wave field (Ap,Am matrices)
            % with a spectrum specified in terms of vertical mode j and
            % frequency $\omega$. This allows us to initialize with a
            % Garrett-Munk spectrum, for example, using code like,
            %
            % ```matlab
            % GM = @(omega,j) E*H(j) .* B(omega);
            % ```
            % 
            % Because the model has limited resolution, there will not
            % necessarily be many modes in a given frequency band. This
            % means that the ensemble may be over a very low number of
            % realization, and thus might not converge to the requested
            % spectrum. For this reason, the option
            % shouldOnlyRandomizeOrientations may be useful. This will only
            % randomize the phases of the waves, while fixing the
            % amplitudes so that the desired spectrum will be achieved.
            %
            % - Topic: Initial conditions — Waves
            % - Declaration: initWavesWithFrequencySpectrum(options)
            % - Parameter ApmSpectrum: function_handle with signature @(omega,j), defaults to a white spectrum.
            % - Parameter shouldOnlyRandomizeOrientations: boolean indicating whether randomness in amplitudes should be eliminated (default 0)
            % - Parameter shouldShowDiagnostics: whether to summarize what just happened (default 0)
            arguments (Input)
                self WVTransform {mustBeNonempty}
                options.ApmSpectrum = @isempty
                options.shouldOnlyRandomizeOrientations (1,1) double {mustBeMember(options.shouldOnlyRandomizeOrientations,[0 1])} = 0
                options.shouldShowDiagnostics (1,1) double {mustBeMember(options.shouldShowDiagnostics,[0 1])} = 0
            end
            self.removeAllWaves;
            optionsArgs = namedargs2cell(options);
            self.addWavesWithFrequencySpectrum(optionsArgs{:});
        end

        function addWavesWithFrequencySpectrum(self,options)
            % add waves with a specified frequency spectrum
            %
            % This allows you to initialize the wave field (Ap,Am matrices)
            % with a spectrum specified in terms of vertical mode j and
            % frequency $\omega$. This allows us to initialize with a
            % Garrett-Munk spectrum, for example, using code like,
            %
            % ```matlab
            % GM = @(omega,j) E*H(j) .* B(omega);
            % ```
            % 
            % Because the model has limited resolution, there will not
            % necessarily be many modes in a given frequency band. This
            % means that the ensemble may be over a very low number of
            % realization, and thus might not converge to the requested
            % spectrum. For this reason, the option
            % shouldOnlyRandomizeOrientations may be useful. This will only
            % randomize the phases of the waves, while fixing the
            % amplitudes so that the desired spectrum will be achieved.
            %
            % - Topic: Initial conditions — Waves
            % - Declaration: addWavesWithFrequencySpectrum(options)
            % - Parameter ApmSpectrum: function_handle with signature @(omega,j), defaults to a white spectrum.
            % - Parameter shouldOnlyRandomizeOrientations: boolean indicating whether randomness in amplitudes should be eliminated (default 0)
            % - Parameter shouldShowDiagnostics: whether to summarize what just happened (default 0)
            arguments (Input)
                self WVTransform {mustBeNonempty}
                options.ApmSpectrum = @isempty
                options.shouldOnlyRandomizeOrientations (1,1) double {mustBeMember(options.shouldOnlyRandomizeOrientations,[0 1])} = 0
                options.shouldShowDiagnostics (1,1) double {mustBeMember(options.shouldShowDiagnostics,[0 1])} = 0
            end

            if isequal(options.ApmSpectrum,@isempty)
                ApmSpectrum = @(omega,j) ones(size(k));
            else
                ApmSpectrum = options.ApmSpectrum;
            end
            
            % We let the spectrum be defined over inertial and wave modes.
            % Perhaps the correct way to do this is to create a flow
            % component with the correct number of modes.
            [Ap_rand,Am_rand,~] = self.waveComponent.randomAmplitudes(shouldOnlyRandomizeOrientations=options.shouldOnlyRandomizeOrientations);
            [Ap_io,Am_io,~] = self.inertialComponent.randomAmplitudes(shouldOnlyRandomizeOrientations=options.shouldOnlyRandomizeOrientations);
            Ap_rand = Ap_rand + Ap_io;
            Am_rand = Am_rand + Am_io;

            % Rather than renormalize the randomized modes, we copy modes
            % once we've touched them. This is good because starting our
            % for-loop at j=1 means we've neglected the j=0 mode, and thus
            % it would be left untouched.
            Ap_ = zeros(self.spectralMatrixSize);
            Am_ = zeros(self.spectralMatrixSize);

            % This loops adds energy to each unique frequency, for each j.
            % It does this by integrating the energy spectrum from the
            % point midway between the frequency to the left, to the
            % point midway between the frequency to the right.
            allIndices = (1:numel(self.Omega)).';
            for j=1:max(self.j) % start at j=1---no barotropic mode!
                indicesForJ = allIndices(self.J == j);
                [sortedOmegas, sortedOmegasIndices] = sort(self.Omega(indicesForJ));
                indicesForJ = indicesForJ(sortedOmegasIndices); % contains the indices back to the Ap matrix
                omegaDiffIndices = find(diff(sortedOmegas) > 0);

                lastIdx = 1;
                omega0 = sortedOmegas(lastIdx); % current omega we will be assigning energy to
                leftDeltaOmega = 0;
                totalEnergyAdded = 0;
                for idx=omegaDiffIndices'
                    indicesForOmega = indicesForJ(lastIdx:idx);
                    nOmegas = length(indicesForOmega);

                    omega1 = sortedOmegas(idx + 1);
                    rightDeltaOmega = (omega1-omega0)/2;

                    energyPerApmComponent = integral(@(omega) ApmSpectrum(omega,j),omega0-leftDeltaOmega,omega0+rightDeltaOmega)/nOmegas/2;
                    Ap_(indicesForOmega) = Ap_rand(indicesForOmega).*sqrt(energyPerApmComponent./(self.Apm_TE_factor(indicesForOmega) ));
                    Am_(indicesForOmega) = Am_rand(indicesForOmega).*sqrt(energyPerApmComponent./(self.Apm_TE_factor(indicesForOmega) ));

                    totalEnergyAdded = totalEnergyAdded + energyPerApmComponent*nOmegas*2;

                    omega0 = omega1;
                    leftDeltaOmega = rightDeltaOmega;
                    lastIdx = idx+1;
                end

                % The right-most omega still needs to be set
                indicesForOmega = indicesForJ(lastIdx:end);
                nOmegas = length(indicesForOmega);

                energyPerApmComponent = integral(@(omega) ApmSpectrum(omega,j),omega0-leftDeltaOmega,omega0)/nOmegas/2;
                Ap_(indicesForOmega) = Ap_rand(indicesForOmega).*sqrt(energyPerApmComponent./(self.Apm_TE_factor(indicesForOmega) ));
                Am_(indicesForOmega) = Am_rand(indicesForOmega).*sqrt(energyPerApmComponent./(self.Apm_TE_factor(indicesForOmega) ));
                totalEnergyAdded = totalEnergyAdded + energyPerApmComponent*nOmegas*2;

                if options.shouldShowDiagnostics == 1
                    totalEnergyDesired = integral( @(omega) ApmSpectrum(omega,j),self.f,sqrt(max(self.N2)));
                    cumulativeEnergy = sum(sum(self.Apm_TE_factor.*(abs(Ap_).^2+abs(Am_).^2)));
                    fprintf('Total energy desired in j=%d is %.5f, total energy added is %.5f. Cumulative is %.5f\n',j,totalEnergyDesired,totalEnergyAdded ,cumulativeEnergy);
                end
            end

            self.throwErrorIfDensityViolation(A0=self.A0,Ap=self.Ap+Ap_,Am=self.Am+Am_,additionalErrorInfo=sprintf('The modes you are setting will cause the fluid state to violate this condition.\n'));
            self.Ap = self.Ap+Ap_;
            self.Am = self.Am+Am_;
        end

        function initWithGMSpectrum(self,options)
            % initialize the wave field following a Garrett-Munk spectrum
            %
            % - Topic: Initial conditions — Waves
            % - Declaration: initWithGMSpectrum(options)
            % - Parameter GMAmplitude: (optional) energy amplitude, default 1
            % - Parameter j_star: (optional) modal roll-off, default 3
            % - Parameter shouldOnlyRandomizeOrientations: (optional) boolean indicating whether randomness in amplitudes should be eliminated (default 0)
            % - Parameter shouldShowDiagnostics: (optional) whether to summarize what just happened (default 0)
            arguments (Input)
                self WVTransform {mustBeNonempty}
                options.GMAmplitude (1,1) double = 1
                options.j_star (1,1) double = 3
                options.shouldOnlyRandomizeOrientations (1,1) double {mustBeMember(options.shouldOnlyRandomizeOrientations,[0 1])} = 0
                options.shouldShowDiagnostics (1,1) double {mustBeMember(options.shouldShowDiagnostics,[0 1])} = 0
            end
            self.removeAllWaves;
            optionsArgs = namedargs2cell(options);
            self.addGMSpectrum(optionsArgs{:});
        end

        function addGMSpectrum(self,options)
            % add waves following a Garrett-Munk spectrum
            %
            % - Topic: Initial conditions — Waves
            % - Declaration: addGMSpectrum(options)
            % - Parameter GMAmplitude: (optional) energy amplitude, default 1
            % - Parameter j_star: (optional) modal roll-off, default 3
            % - Parameter shouldOnlyRandomizeOrientations: (optional) boolean indicating whether randomness in amplitudes should be eliminated (default 0)
            % - Parameter shouldShowDiagnostics: (optional) whether to summarize what just happened (default 0)
            arguments (Input)
                self WVTransform {mustBeNonempty}
                options.GMAmplitude (1,1) double = 1
                options.j_star (1,1) double = 3
                options.shouldOnlyRandomizeOrientations (1,1) double {mustBeMember(options.shouldOnlyRandomizeOrientations,[0 1])} = 0
                options.shouldShowDiagnostics (1,1) double {mustBeMember(options.shouldShowDiagnostics,[0 1])} = 0
            end
            L_gm = 1.3e3; % thermocline exponential scale, meters
            invT_gm = 5.2e-3; % reference buoyancy frequency, radians/seconds
            E_gm = 6.3e-5; % non-dimensional energy parameter
            E = options.GMAmplitude*L_gm*L_gm*L_gm*invT_gm*invT_gm*E_gm;

            j_star = options.j_star;
            Nmax = sqrt(max(self.N2));

            % Compute the proper vertical function normalization
            H = (j_star+(1:1024)).^(-5/2);
            H_norm = 1/sum(H);

            % Do the same for the frequency function.
            B_norm = 1/atan(sqrt(Nmax*Nmax/(self.f*self.f)-1));

            H = @(j) H_norm*((j+j_star).^(-5/2));
            B = @(omega) B_norm*self.f./(omega.*sqrt(omega.*omega-self.f*self.f));
            GM = @(omega,j) E*H(j) .* B(omega);

            self.addWavesWithFrequencySpectrum(ApmSpectrum=GM,shouldOnlyRandomizeOrientations=options.shouldOnlyRandomizeOrientations,shouldShowDiagnostics=0);

            if options.shouldShowDiagnostics == 1
                igwEnergy = self.Apm_TE_factor.*(abs(self.Ap).^2+abs(self.Am).^2);
                totalEnergyActual = sum(sum(igwEnergy));
                fprintf('After distributing energy across frequency and mode, you have %.2f%% of reference GM energy.\n',100*totalEnergyActual/E);
                if options.shouldOnlyRandomizeOrientations == 0
                    fprintf('Because you are allowing amplitudes to be randomized, primarly due to the amplitudes of the most energetic modes.\n')
                end


                verticalModeTotalEnergy = 0;
                for j=1:max(self.j)
                    verticalModeTotalEnergy = verticalModeTotalEnergy + integral( @(omega) GM(omega,j),self.f,Nmax);
                end
                missingVerticalModeEnergy = E - verticalModeTotalEnergy;
                fprintf('\tYou are missing %.2f%% of the energy due to limited vertical modes.\n',100*missingVerticalModeEnergy/E);
                fprintf('\tYou are missing %.2f%% of the energy due to limited horizontal resolution.\n',100*(E-totalEnergyActual-missingVerticalModeEnergy)/E);

                
                fprintf('Due to restricted domain size:\n');
                fprintf('\tThe inertial oscillations contain %.2f%% of the total energy\n',100*sum(sum(igwEnergy(self.Kh == 0)))/totalEnergyActual);
                fprintf('\tThe inertial oscillation j=1 mode contains %.2f%% of the total energy\n',100*igwEnergy(self.J==1 & self.Kh == 0)/totalEnergyActual);
            end

        end

        function initWithAlternativeSpectrum(self,options)
            % initialize with an alternative formulation of the GM spectrum in the wavenumber domain.
            %
            % This only initializes the wave components, A0 is left untouched.
            %
            % - Topic: Initial conditions — Waves
            % - Declaration: initWithAlternativeSpectrum(options)
            % - Parameter GMAmplitude: (optional) default 1
            % - Parameter j_star: (optional) default 3
            % - Parameter slope: (optional) default 1
            arguments (Input)
                self WVTransform {mustBeNonempty}
                options.GMAmplitude (1,1) double =1
                options.j_star (1,1) double = 3
                options.slope (1,1) double = 1
                %options.shouldRandomizeAmplitude = 1
            end

            j_star=options.j_star;
            slope=options.slope;
            GMAmplitude=options.GMAmplitude;

            % GM Parameters
            L_gm = 1.3e3; % thermocline exponential scale, meters
            invT_gm = 5.2e-3; % reference buoyancy frequency, radians/seconds
            E_gm = 6.3e-5; % non-dimensional energy parameter
            E_T = L_gm*L_gm*L_gm*invT_gm*invT_gm*E_gm*GMAmplitude;

            % Compute the proper vertical function normalization
            M = @(j) (j_star.^2 +(j).^2).^((-5/4));
            M_norm = sum(M(1:1024));
            M= @(j) ((j_star.^2 +(j).^2).^((-5/4)))/M_norm;

            % sanity check to confirm this is 1
            % sum(M(1:1024))

            % Definir a função anônima fun (k?)
            B = @(k, j) (1./(k.^2.* self.Lr2(j+1) + 1).^(1 * slope)).*sqrt(self.Lr2(j+1));

            % Definir a função anônima B_norm que integra fun em relação a k
            B_norm = ones(self.Nj,1);
            for jind=(2:self.Nj)                
                B_norm(jind) = integral(@(k) B(k, self.j(jind)), 0, 1);
            end

            B = @(k, j) (1./(k.^2.* self.Lr2(j+1) + 1).^(1 * slope)).*sqrt(self.Lr2(j+1))/B_norm(j+1);

            % Sanity check to confirm that the integrals are now normalized
            % for jind=(2:self.Nj)                
            %     integral(@(k) B(k, self.j(jind)), 0, 1)
            % end
      
            % Definir a função model_spectrum
            model_spectrum = @(k, j) (E_T) * B(k, j) * M(j);

            [self.Ap,self.Am,~] = self.waveComponent.randomAmplitudesWithSpectrum(ApmSpectrum= model_spectrum,shouldOnlyRandomizeOrientations=1);
        end

    end

    methods (Hidden=true)
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

            [A_plus,phi_plus,linearIndex] = WVTransform.extractNonzeroWaveProperties(A_p);
            omega_plus = self.Omega(linearIndex);
            mode_plus = J(linearIndex);
            alpha_plus = atan2(L(linearIndex),K(linearIndex));
            k_plus = K(linearIndex);
            l_plus = L(linearIndex);

            [A_minus,phi_minus,linearIndex] = WVTransform.extractNonzeroWaveProperties(A_m);
            omega_minus = -self.Omega(linearIndex);
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
    end

    methods (Static, Hidden=true)
        function propertyAnnotations = propertyAnnotationsForInternalGravityWaveComponent(options)
            % return array of WVVariableAnnotation instances initialized by default
            %
            % This function creates annotations for the built-in variables supported by
            % the WVTransform.
            %
            % - Topic: Internal
            % - Declaration: operations = defaultVariableAnnotations()
            % - Returns operations: array of WVVariableAnnotation instances
            arguments
                options.spectralDimensionNames = {'j','kl'}
            end
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);

            annotation = WVVariableAnnotation('waveEnergy',{},'m^3 s^{-2}', 'total energy, waves');
            annotation.isVariableWithLinearTimeStep = 0;
            annotation.isVariableWithNonlinearTimeStep = 1;
            propertyAnnotations(end+1) = annotation;

            propertyAnnotations(end+1) = WVPropertyAnnotation('Omega',{'k','l','j'},'rad s^{-1}', 'frequency of oscillation of the linear waves', detailedDescription='- topic: Domain Attributes');


            propertyAnnotations(end+1) = CANumericProperty('ApU',options.spectralDimensionNames,'', 'matrix component that multiplies $$\tilde{u}$$ to compute $$A_p$$.',isComplex=1);
            propertyAnnotations(end+1) = CANumericProperty('ApV',options.spectralDimensionNames,'', 'matrix component that multiplies $$\tilde{v}$$ to compute $$A_p$$.',isComplex=1);
            propertyAnnotations(end+1) = CANumericProperty('ApN',options.spectralDimensionNames,'s^{-1}', 'matrix component that multiplies $$\tilde{\eta}$$ to compute $$A_p$$.',isComplex=0);
            propertyAnnotations(end+1) = CANumericProperty('AmU',options.spectralDimensionNames,'', 'matrix component that multiplies $$\tilde{u}$$ to compute $$A_m$$.',isComplex=1);
            propertyAnnotations(end+1) = CANumericProperty('AmV',options.spectralDimensionNames,'', 'matrix component that multiplies $$\tilde{v}$$ to compute $$A_m$$.',isComplex=1);
            propertyAnnotations(end+1) = CANumericProperty('AmN',options.spectralDimensionNames,'s^{-1}', 'matrix component that multiplies $$\tilde{\eta}$$ to compute $$A_m$$.',isComplex=0);


            propertyAnnotations(end+1) = CANumericProperty('UAp',options.spectralDimensionNames,'', 'matrix component that multiplies $$A_p$$ to compute $$\tilde{u}$$.',isComplex=1);
            propertyAnnotations(end+1) = CANumericProperty('UAm',options.spectralDimensionNames,'', 'matrix component that multiplies $$A_m$$ to compute $$\tilde{u}$$.',isComplex=1);
            propertyAnnotations(end+1) = CANumericProperty('VAp',options.spectralDimensionNames,'', 'matrix component that multiplies $$A_p$$ to compute $$\tilde{v}$$.',isComplex=1);
            propertyAnnotations(end+1) = CANumericProperty('VAm',options.spectralDimensionNames,'', 'matrix component that multiplies $$A_m$$ to compute $$\tilde{v}$$.',isComplex=1);
            propertyAnnotations(end+1) = CANumericProperty('WAp',options.spectralDimensionNames,'', 'matrix component that multiplies $$A_p$$ to compute $$\tilde{w}$$.',isComplex=1);
            propertyAnnotations(end+1) = CANumericProperty('WAm',options.spectralDimensionNames,'', 'matrix component that multiplies $$A_m$$ to compute $$\tilde{w}$$.',isComplex=1);
            propertyAnnotations(end+1) = CANumericProperty('NAp',options.spectralDimensionNames,'s', 'matrix component that multiplies $$A_p$$ to compute $$\tilde{\eta}$$.',isComplex=0);
            propertyAnnotations(end+1) = CANumericProperty('NAm',options.spectralDimensionNames,'s', 'matrix component that multiplies $$A_m$$ to compute $$\tilde{\eta}$$.',isComplex=0);

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