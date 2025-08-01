%% Load the file

file = 'ForcedDissipativeQG-particles-512.nc';

ncfile = NetCDFFile(file);
tDim = ncfile.readVariables('t');
k_f = ncfile.attributes('k_f');
k_r = ncfile.attributes('k_r');
L_damp = ncfile.readVariables('L_damp');

iTime = 350;
[F0_psi,F_psi,qgpv] = ncfile.readVariablesAtIndexAlongDimension('t',iTime,'F0_psi','F_psi','qgpv');
eta = mean(mean(qgpv .* F_psi)); % 1/s * 1/s^2 = 1/s^3
fprintf('Enstrophy forcing, computed spatially: %g 1/s^3\n',eta);

wvt = WVTransform.waveVortexTransformFromFile(file,iTime=iTime);

% Muliply A0 by this factor and you get QGPV
PVFactor = wvt.A0_QGPV_factor; % real, units of 1/(m s)
FZ_forcing = PVFactor.*F0_psi;
FZ_damping = PVFactor.*(L_damp.*wvt.A0);

enstrophyInertial = -wvt.enstrophyFlux();
enstrophyForcing = -PVFactor .* real( FZ_forcing .* conj(wvt.A0) );
enstrophyDamping = PVFactor .* real( FZ_damping .* conj(wvt.A0) );


[FZ_inertial,FZ_forcing,FZ_damping] = wvt.transformToRadialWavenumber(enstrophyInertial,enstrophyForcing,enstrophyDamping);

figure
plot(wvt.kRadial,cumsum(FZ_inertial),LineWidth=2), xlog, hold on
plot(wvt.kRadial,cumsum(FZ_forcing),LineWidth=2)
plot(wvt.kRadial,cumsum(FZ_damping),LineWidth=2)
plot(wvt.kRadial,cumsum(FZ_forcing+FZ_damping),LineWidth=2)
legend('inertial','forcing','damping','forcing+damping')
% vlines([k_f k_r],'g--')

return;

ncfile = NetCDFFile(file);
tDim = ncfile.readVariables('t');
k_f = ncfile.attributes('k_f');
k_r = ncfile.attributes('k_r');
L_damp = ncfile.readVariables('L_damp');

% Muliply A0 by this factor and you get QGPV
PVFactor = -wvt.Omega .* wvt.Omega / (wvt.h * wvt.f); % real, units of 1/(m s)
TotalEnstrophy = (PVFactor .* wvt.A0) .* conj(PVFactor.*wvt.A0);
Z0_radial = wvt.transformToRadialWavenumber(TotalEnstrophy);

figure, plot(wvt.kRadial,Z0_radial), xlog, ylog, hold on
ylabel('1/s^2')
xlabel('1/m')
title('enstrophy spectral density')
vlines([k_f,k_r],'g--')


timeIndices = 1:1;

TE = zeros(length(timeIndices),length(Z0_radial));
ZF_inertial = zeros(size(TE));
ZF_forcing = zeros(size(TE));
ZF_damping = zeros(size(TE));

for iT = 1:length(timeIndices)
    iTime = timeIndices(iT);
    wvt = WVTransform.waveVortexTransformFromFile(file,iTime=iTime);

    F0_psi = ncfile.readVariablesAtIndexAlongDimension('t',iTime,'F0_psi');

    TotalEnstrophy = (PVFactor.*PVFactor) .* (wvt.A0.*conj(wvt.A0)); % 1/s^2
    EnstrophyFlux = wvt.enstrophyFlux(); % 1/s^3
    EnstrophyFluxForcing = 2*(PVFactor.*PVFactor) .* real( F0_psi .* conj(wvt.A0) ); % 1/(m s) * m/s * 1/s = 1/s^3
    EnstrophyFluxDamping = 2*PVFactor .* real( (L_damp.*(PVFactor.*wvt.A0)) .* conj(PVFactor.*wvt.A0) );

    [TotalEnstrophyRadial,EnstrophyFluxRadial,EnstrophyFluxForcingRadial,EnstrophyFluxDampingRadial] = wvt.transformToRadialWavenumber(TotalEnstrophy,EnstrophyFlux,EnstrophyFluxForcing,EnstrophyFluxDamping);

    TE(iT,:) = TotalEnstrophyRadial;
    ZF_inertial(iT,:) = EnstrophyFluxRadial;
    ZF_forcing(iT,:) = EnstrophyFluxForcingRadial;
    ZF_damping(iT,:) = EnstrophyFluxDampingRadial;
end

dk = (wvt.kRadial(2)-wvt.kRadial(1));

figure
plot(wvt.kRadial,mean(TE,1)/dk), xlog, ylog, hold on
ylabel('m^3/s^2')
xlabel('1/m')
title('enstrophy spectrum, spectral density')
vlines([k_f,k_r],'g--')

figure
% plot(wvt.kRadial,mean(-ZF_inertial,1)), xlog, hold on
% plot(wvt.kRadial,mean(ZF_forcing,1))
plot(wvt.kRadial,mean(ZF_damping,1))
title('enstrophy flux into a given wavenumber')
ylabel('1/s^3')
legend('inertial','forcing','damping')

% The flux rate into a given wavenumber is just the amount removed by
% damping (of course) in steady-state. The total integral of which, is the
% forcing at the forcing wavenumber.
return
timeIndices = 1:length(tDim);
EF_forcingT = zeros(length(timeIndices),1);
for iT = 1:length(timeIndices)
    iTime = timeIndices(iT);
    wvt = WVTransform.waveVortexTransformFromFile(file,iTime=iTime);
    F0_psi = ncfile.readVariablesAtIndexAlongDimension('t',iTime,'F0_psi');
    EnstrophyFluxForcing = 2*wvt.A0_TE_factor.*real( F0_psi .* conj(wvt.A0) ); % m/s^2 * m/s * m = m^3/s^3
    EF_forcingT(iT) = -sum(EnstrophyFluxForcing(:));
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

%%
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

%%
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
