%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions and initialize a WVTransform.
% The 'h' parameter is the equivalent depth, and 0.80 m is a typical value
% for the first baroclinic mode.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lxy = 50e3;
Nxy = 256;
latitude = 25;

wvt = WVTransformSingleMode([Lxy, Lxy], [Nxy, Nxy], h=0.8, latitude=latitude);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize a non-linear flux operator---in this case, the NonlinearFlux
% subclass is specific to QG forced-dissipative turbulence.
%
% We set a damping wavenumber, k_r, a forcing wavenumber, k_f, and a target
% rms velocity, u_rms. Both k_r and u_rms are *approximate*.
%
% initialPV can be set to none (which is useful when restarting),
% 'narrow-band', or 'full-spectrum'. narrow-band will take a while to spin
% up, but conversely, full-spectrum will take a while to adjust.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dk = (wvt.k(2)-wvt.k(1));
k_f = 30*dk;
k_r = 4*dk;
u_rms = 0.05;

fdFlux = SingleModeForcedDissipativeQGPVEMasked(wvt,k_f=k_f,k_r=k_r,u_rms=u_rms,initialPV='full-spectrum');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Once a WVTransform and a NonlinearFlux operator have been
% initialized, we can now initialize a model.
% 
% Here we choose a relatively small time step, as we expect energy to build
% up over time.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model = WVModel(wvt,nonlinearFlux=fdFlux);
model.setupIntegrator(deltaT=0.5*model.nonlinearFluxOperation.dampingTimeScale,outputInterval=86400);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Set an output file, set the variables that we want written to file, and
% integrate. How long to integrate for? We need to wait until energy
% reaches steady-state. You can always keep integrating to be certain.
%
% For these particular settings, I find 150 days is a good initial
% integration.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% model.createNetCDFFileForModelOutput(sprintf('ForcedDissipativeQG-spinup-%d.nc',Nxy),shouldOverwriteExisting=1);
% model.setNetCDFOutputVariables('A0','psi','zeta_z','F_psi','F0_psi');
model.integrateToTime(50*86400);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
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
plot(wvt.kRadial,fdFlux.model_spectrum(wvt.kRadial))
ylabel('m^3/s^2')
xlabel('1/m')
title('horizontal velocity spectrum')
vlines([k_f,k_r],'g--')