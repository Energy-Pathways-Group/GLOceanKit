file = 'ForcedDissipativeQG-256-restart.nc';
iTime = 10;
wvt = WaveVortexTransform.transformFromFile(file,iTime=Inf);
TotalEnergy = (wvt.A0_TE_factor/wvt.h) .* (wvt.A0.*conj(wvt.A0));
E_damp_radial = wvt.transformToRadialWavenumber(TotalEnergy);

timeIndices = 1:10;

TE = zeros(length(timeIndices),length(E_damp_radial));
EF_inertial = zeros(size(TE));
EF_forcing = zeros(size(TE));
EF_damping = zeros(size(TE));

ncfile = NetCDFFile(file);
tDim = ncfile.readVariables('t');
L_damp = ncfile.readVariables('L_damp');
k_f = ncfile.attributes('k_f');
k_r = ncfile.attributes('k_r');

for iT = 1:length(timeIndices)
    iTime = timeIndices(iT);
    wvt = WaveVortexTransform.transformFromFile(file,iTime=iTime);

    F0_psi = ncfile.readVariablesAtIndexAlongDimension('t',iTime,'F0_psi');

    TotalEnergy = (wvt.A0_TE_factor/wvt.h) .* (wvt.A0.*conj(wvt.A0)); % m^2/s^3

    [~,~,EnergyFluxInertial] = wvt.energyFlux(); % m^3/s^3
    EnergyFluxForcing = 2*wvt.A0_TE_factor.*real( F0_psi .* conj(wvt.A0) ); % m/s^2 * m/s * m = m^3/s^3
    EnergyFluxDamping = 2*wvt.A0_TE_factor.*real( (L_damp.*wvt.A0) .* conj(wvt.A0) ); % m/s^2 * m/s * m = m^3/s^3

    [TotalEnergyRadial,EnergyFluxInertialRadial,EnergyFluxForcingRadial,EnergyFluxDampingRadial] = wvt.transformToRadialWavenumber(TotalEnergy,EnergyFluxInertial,EnergyFluxForcing,EnergyFluxDamping);

    TE(iT,:) = TotalEnergyRadial;
    EF_inertial(iT,:) = EnergyFluxInertialRadial;
    EF_forcing(iT,:) = EnergyFluxForcingRadial;
    EF_damping(iT,:) = EnergyFluxDampingRadial;
end

dk = (wvt.kRadial(2)-wvt.kRadial(1));

figure
plot(wvt.kRadial,mean(TE,1)/dk), xlog, ylog, hold on
ylabel('m^3/s^2')
xlabel('1/m')
title('horizontal velocity spectrum, spectral density')
vlines([k_f,k_r],'g--')

figure
plot(wvt.kRadial,mean(-EF_inertial,1)/wvt.h), xlog, hold on
plot(wvt.kRadial,mean(EF_forcing,1)/wvt.h)
plot(wvt.kRadial,mean(EF_damping,1)/wvt.h)
title('energy flux into a given wavenumber')
ylabel('m^2/s^3')

% The flux rate into a given wavenumber is just the amount removed by
% damping (of course) in steady-state. The total integral of which, is the
% forcing at the forcing wavenumber.

timeIndices = 1:length(tDim);
EF_forcingT = zeros(length(timeIndices),1);
for iT = 1:length(timeIndices)
    iTime = timeIndices(iT);
    wvt = WaveVortexTransform.transformFromFile(file,iTime=iTime);
    F0_psi = ncfile.readVariablesAtIndexAlongDimension('t',iTime,'F0_psi');
    EnergyFluxForcing = 2*wvt.A0_TE_factor.*real( F0_psi .* conj(wvt.A0) ); % m/s^2 * m/s * m = m^3/s^3
    EF_forcingT(iT) = -sum(EnergyFluxForcing(:));
end
figure
plot(tDim/86400,EF_forcingT/wvt.h)
ylabel('energy forcing (m^2/s^3)')
xlabel('time (days)')

r = ncfile.attributes('r');
figure
plot(tDim/86400, sqrt(EF_forcingT/(r*wvt.h*(2*pi)^2)))
ylabel('energy forcing (m/s)')
xlabel('time (days)')

return

    dk = (wvt.kRadial(2)-wvt.kRadial(1));


figure
plot(wvt.kRadial,EkT/dk), xlog, ylog, hold on
ylabel('m^3/s^2')
xlabel('1/m')
title('horizontal velocity spectrum, spectral density')
vlines([k_f,k_r],'g--')


[K,L,~] = ndgrid(wvt.k,wvt.l,wvt.j);
Kh = sqrt(K.*K + L.*L);
kRadial = wvt.kRadial;
AA = ~(wvt.maskForAliasedModes(jFraction=1));
energyFlux = zeros(length(kRadial),1);
for iK=1:length(kRadial)
    A0Mask = AA;
    A0Mask(Kh > kRadial(iK)-dk/2 & Kh < kRadial(iK)+dk/2) = 0;
    [Ep,Em,E0] = wvt.energyFluxWithMasks(zeros(size(wvt.A0)),zeros(size(wvt.A0)),A0Mask);
    energyFlux(iK) = sum(E0(:));
end
figure, plot(kRadial,energyFlux/wvt.h), xlog
title('energy flux into a given wavenumber')
ylabel('m^2/s^3')


ncfile = NetCDFFile(file);
tDim = ncfile.readVariables('t');
[psi,F_psi,F0_psi,zeta_z,A0,L_damp] = ncfile.readVariablesAtIndexAlongDimension('t',iTime,'psi','F_psi','F0_psi','zeta_z','A0','L_damp');
epsilon = mean(mean(psi .* F_psi)) % m^2/s * 1/s^2 = m^2/s^3
eta = mean(mean(zeta_z .* F_psi)) % 1/s * 1/s^2 = 1/s^3

E_force = 2*wvt.A0_TE_factor.*real( F0_psi .* conj(A0) ); % m/s^2 * m/s * m = m^3/s^3
E_force_radial = wvt.transformToRadialWavenumber(E_force);
sum(E_force(:))/wvt.h

E_damp = 2*wvt.A0_TE_factor.*real( (L_damp.*A0) .* conj(A0) ); % m/s^2 * m/s * m = m^3/s^3
E_damp_radial = wvt.transformToRadialWavenumber(E_damp);
sum(E_damp(:))/wvt.h


figure
plot(kRadial,energyFlux/wvt.h), hold on
plot(kRadial,-E_damp_radial/wvt.h), xlog


figure, plot(kRadial,cumsum(energyFlux)/wvt.h), xlog
hold on
plot(kRadial,cumsum(-E_damp_radial)/wvt.h)
plot(kRadial,cumsum(-E_force_radial)/wvt.h)
