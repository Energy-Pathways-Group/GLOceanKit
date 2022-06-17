model = WaveVortexModel.modelFromFile('ForcedDissipativeQG.nc');

deltaT = (wvt.x(2)-wvt.x(1))*0.25/(pi*model.nonlinearFlux.u_rms);
model.setupIntegrator(deltaT=deltaT, outputInterval=86400);
model.createNetCDFFileForModelOutput('ForcedDissipativeQG-restart.nc',shouldOverwriteExisting=1);
model.integrateToTime(450*86400);

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
vlines([k_f,k_r],'g--')