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
    
    [K,L,~] = ndgrid(self.k,self.l,self.j);
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
        A_plus = A.*WaveVortexModel.generateHermitianRandomMatrix( size(K) );
        A_minus = A.*WaveVortexModel.generateHermitianRandomMatrix( size(K) );

        self.offgridModes.U_ext = sqrt(2*GM3Dext./self.offgridModes.h_ext).*randn( size(self.offgridModes.h_ext) );
        self.offgridModes.PrecomputeExternalWaveCoefficients();                
    else
        % Randomize phases, but keep unit length
        A_plus = WaveVortexModel.generateHermitianRandomMatrix( size(K), excludeNyquist );
        A_minus = WaveVortexModel.generateHermitianRandomMatrix( size(K), excludeNyquist );

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