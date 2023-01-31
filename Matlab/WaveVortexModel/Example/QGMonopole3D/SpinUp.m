%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Specify the problem dimensions and initialize a WVTransform.
% The 'h' parameter is the equivalent depth, and 0.80 m is a typical value
% for the first baroclinic mode.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lxy = 50e3;
Nxy = 128;
Nz = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lz = 4000;
N0 = 3*2*pi/3600; % buoyancy frequency at the surface, radians/seconds
L_gm = 1300; % thermocline exponential scale, meters
N2 = @(z) N0*N0*exp(2*z/L_gm);
wvt = WVTransformHydrostatic([Lxy, Lxy, Lz], [Nxy, Nxy, Nz], N2=N2,latitude=25);

outputVar = WVVariableAnnotation('zeta_z',{'x','y','z'},'1/s^2', 'vertical component of relative vorticity');
f = @(wvt) wvt.transformToSpatialDomainWithF(-(wvt.g/wvt.f) * ((wvt.K).^2 +(wvt.L).^2) .* wvt.A0t);
wvt.addOperation(WVOperation('zeta_z',outputVar,f));

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
k_f = 15*dk;
k_r = 4*dk;
u_rms = 0.05;

fdFlux = ForcedDissipativeQGPVE(wvt,k_f=k_f,k_r=k_r,u_rms=u_rms,initialPV='narrow-band');
% fdFlux = ForcedDissipativeQGPVE(wvt,k_f=k_f,k_r=k_r,u_rms=u_rms,initialPV='full-spectrum');

ssu = wvt.seaSurfaceU;
ssv = wvt.seaSurfaceV;
ssh = wvt.seaSurfaceHeight;

mean(mean(ssu.^2 + ssv.^2))
wvt.g*mean(mean(ssh.^2))/wvt.h(2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
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
%%
% Set an output file, set the variables that we want written to file, and
% integrate. How to integrate for? We need to wait until energy reaches
% steady-state. You can always keep integrating to be certain.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model.createNetCDFFileForModelOutput(sprintf('ForcedDissipativeQG-spinup-%d.nc',Nxy),shouldOverwriteExisting=1);
model.setNetCDFOutputVariables('A0','psi','zeta_z','F_psi','F0_psi');
model.integrateToTime(25*86400);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% At the end of the integration, let's make a plot showing the relative
% vorticity and resulting energy spectrum.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EkT = wvt.transformToRadialWavenumber(wvt.A0_TE_factor.* (wvt.A0.*conj(wvt.A0)));

zeta_z = wvt.zeta_z;
figure(Position=[100 100 1000 1000])
tiledlayout('flow')

nexttile
pcolor(wvt.x/1000,wvt.y/1000,zeta_z(:,:,end)), shading interp
colormap("gray")
xlabel('km'), ylabel('km')

nexttile
pcolor(wvt.y/1000,wvt.z/1000,squeeze(zeta_z(1,:,:)).'), shading interp
colormap("gray")
xlabel('km'), ylabel('km')

nexttile
plot(wvt.kRadial,EkT(:,1)/(wvt.kRadial(2)-wvt.kRadial(1))), xlog, ylog, hold on
plot(wvt.kRadial,EkT(:,2)/(wvt.kRadial(2)-wvt.kRadial(1)))
plot(wvt.kRadial,EkT(:,3)/(wvt.kRadial(2)-wvt.kRadial(1)))
plot(wvt.kRadial,fdFlux.model_spectrum(wvt.kRadial))
ylabel('m^3/s^2')
xlabel('1/m')
title('horizontal velocity spectrum')
vlines([k_f,k_r],'g--')
legend('mode 0','mode 1', 'mode 2', 'model spectrum', 'k_r', 'k_f')

nexttile
plot(wvt.j, sum(EkT,1) )

