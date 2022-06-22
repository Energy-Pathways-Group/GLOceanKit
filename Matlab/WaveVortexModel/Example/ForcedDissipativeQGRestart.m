model = WaveVortexModel.modelFromFile('ForcedDissipativeQG-256.nc',shouldDoubleResolution=1,restartIndex=150);
wvt = model.wvt;

deltaT = (wvt.x(2)-wvt.x(1))*0.25/(pi*model.nonlinearFlux.u_rms);
model.setupIntegrator(deltaT=deltaT, outputInterval=86400);
model.createNetCDFFileForModelOutput(sprintf('ForcedDissipativeQG-%d-restart.nc',Nxy),shouldOverwriteExisting=1);
model.setNetCDFOutputVariables('A0','psi','zeta_z','F_psi','F0_psi');
model.integrateToTime(wvt.t + 10*86400);

EkT = wvt.transformToRadialWavenumber((wvt.A0_TE_factor/wvt.h) .* (wvt.A0.*conj(wvt.A0)));

figure
plot(wvt.kRadial,EkT/(wvt.kRadial(2)-wvt.kRadial(1))), xlog, ylog, hold on
ylabel('m^3/s^2')
xlabel('1/m')
title('initial horizontal velocity spectrum (randomized phases)')
vlines([k_f,k_r],'g--')

return;

model = WaveVortexModel.modelFromFile('ForcedDissipativeQG-restart.nc',shouldDoubleResolution=1);
wvt = model.wvt;
deltaT = (wvt.x(2)-wvt.x(1))*0.25/(pi*model.nonlinearFlux.u_rms);
model.setupIntegrator(deltaT=deltaT, outputInterval=86400);
model.integrateToTime(wvt.t + 5*86400);

model = WaveVortexModel(wvt,nonlinearFlux=model.nonlinearFlux);
model.setupIntegrator(deltaT=deltaT, outputInterval=86400);
model.createNetCDFFileForModelOutput('ForcedDissipativeQG-high-res.nc',shouldOverwriteExisting=1);
model.integrateToTime(wvt.t + 5*86400);

u2 = wvt.u.^2 + wvt.v.^2;
u_max = max(sqrt(u2(:)))
u_rms = sqrt(mean(u2(:)))

% figure, pcolor(wvt.x,wvt.y,wvt.ssh.'), shading interp

% model.integrateToNextOutputTime();
A02 = (wvt.A0_TE_factor/wvt.h) .* (wvt.A0.*conj(wvt.A0));
u_rms_alt = sqrt(2*sum(A02(:)))
GeostrophicEnergyK = wvt.transformToRadialWavenumber(A02);
dk = wvt.kRadial(2)-wvt.kRadial(1);
figure
plot(wvt.kRadial,GeostrophicEnergyK/dk), xlog, ylog
ylabel('m^3/s^2')
xlabel('1/m')
vlines([model.nonlinearFlux.k_f,model.nonlinearFlux.k_r],'g--')

[K,L,~] = ndgrid(wvt.k,wvt.l,wvt.j);
Kh = sqrt(K.*K + L.*L);
kRadial = wvt.kRadial;
AA = ~(wvt.MaskForAliasedModes(jFraction=1));
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

dy = log(GeostrophicEnergyK(6))-log(GeostrophicEnergyK(2));
dx = log(kRadial(6))-log(kRadial(2));
dy/dx