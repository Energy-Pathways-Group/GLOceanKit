%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lxy = 50e3;
Nxy = 256;
latitude = 25;

wvt = WaveVortexTransformSingleMode([Lxy, Lxy], [Nxy, Nxy], h=0.8, latitude=latitude);

% Forcing parameters
dk = (wvt.k(2)-wvt.k(1));
k_f = 30*dk;
k_r = 4*dk;
u_rms = 0.05;


fdFlux = SingleModeForcedDissipativeQGPVEMasked(wvt,k_f=k_f,k_r=k_r,u_rms=u_rms,initialPV='narrow-band');
% Record the initial spectrum that was generated
Ek0 = wvt.transformToRadialWavenumber((wvt.A0_TE_factor/wvt.h) .* (wvt.A0.*conj(wvt.A0)));

model = WaveVortexModel(wvt,nonlinearFlux=fdFlux);
model.setupIntegrator(deltaT=0.1*(wvt.x(2)-wvt.x(1))/u_rms,outputInterval=86400);
model.createNetCDFFileForModelOutput(sprintf('ForcedDissipativeQG-%d.nc',Nxy),shouldOverwriteExisting=1);
model.setNetCDFOutputVariables('A0','psi','zeta_z','F_psi','F0_psi');
model.integrateToTime(60*86400);

EkT = wvt.transformToRadialWavenumber((wvt.A0_TE_factor/wvt.h) .* (wvt.A0.*conj(wvt.A0)));

figure
plot(wvt.kRadial,Ek0/(wvt.kRadial(2)-wvt.kRadial(1))), xlog, ylog, hold on
plot(wvt.kRadial,EkT/(wvt.kRadial(2)-wvt.kRadial(1)))
plot(wvt.kRadial,fdFlux.model_spectrum(wvt.kRadial))
ylabel('m^3/s^2')
xlabel('1/m')
title('initial horizontal velocity spectrum (randomized phases)')
vlines([k_f,k_r],'g--')

return

u2 = wvt.u.^2 + wvt.v.^2;
u_max = max(sqrt(u2(:)))
u_rms = sqrt(mean(u2(:)))

% figure, pcolor(log10(abs(fftshift(fdFlux.damp)))), shading interp
% figure, pcolor(log10(abs(fftshift(fdFlux.F)))), shading interp
% 
% DampK = wvt.transformToRadialWavenumber(fdFlux.damp);