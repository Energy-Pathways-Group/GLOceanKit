%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize a new WaveVortexModel (which includes a WaveVortexTransform
% and a NonlinearFluxOperation) from existing output. We will start from
% the final time-point, and double the resolution.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model = WaveVortexModel.modelFromFile('ForcedDissipativeQG-spinup-256.nc',shouldDoubleResolution=1,restartIndex=Inf);
wvt = model.wvt;

model.setupIntegrator(deltaT=0.5*model.nonlinearFlux.dampingTimeScale,outputInterval=86400);
model.createNetCDFFileForModelOutput('ForcedDissipativeQG-spinup-512.nc',shouldOverwriteExisting=1);
model.setNetCDFOutputVariables('A0','psi','zeta_z','F_psi','F0_psi');
model.integrateToTime(wvt.t + 10*86400);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% At the end of the integration, let's make a plot showing the relative
% vorticity and resulting energy spectrum.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EkT = wvt.transformToRadialWavenumber((wvt.A0_TE_factor/wvt.h) .* (wvt.A0.*conj(wvt.A0)));

figure(Position=[100 100 1000 400])
subplot(1,2,1)
pcolor(wvt.X/1000,wvt.Y/1000,wvt.zeta_z), shading interp
colormap("gray")
xlabel('km'), ylabel('km')
subplot(1,2,2)
plot(wvt.kRadial,EkT/(wvt.kRadial(2)-wvt.kRadial(1))), xlog, ylog, hold on
ylabel('m^3/s^2')
xlabel('1/m')
title('horizontal velocity spectrum')
k_f = model.ncfile.attributes('k_f');
k_r = model.ncfile.attributes('k_r');
vlines([k_f,k_r],'g--')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% What is the enstrophy cascade rate time scale? This is useful for setting
% the output time when we add particles.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tDim = model.ncfile.readVariables('t');
[psi,F_psi] = model.ncfile.readVariablesAtIndexAlongDimension('t',length(tDim),'psi','F_psi');
eta = mean(mean(wvt.zeta_z .* wvt.F_psi)); % 1/s * 1/s^2 = 1/s^3
fprintf('Enstrophy cascade rate has a time scale of %f days.\n',(eta^(-1/3))/86400);