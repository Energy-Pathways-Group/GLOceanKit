%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Specify the problem dimensions and initialize a WVTransform.
% The 'h' parameter is the equivalent depth, and 0.80 m is a typical value
% for the first baroclinic mode.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lxy = 50e3;
Nxy = 512;
latitude = 25;

wvt = WVTransformSingleMode([Lxy, Lxy], [Nxy, Nxy], h=0.8, latitude=latitude);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
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

k_f = 30*wvt.dk;
k_r = 4*wvt.dk;
u_rms = 0.05;
uv_damp = (pi^2)*u_rms; % Damping needs to know the maximum value of u; used to set the small scale damping

fdFlux = WVNonlinearFluxQGForced(wvt,uv_damp=uv_damp);
model_spectrum2D = fdFlux.setNarrowBandForcing(initialPV='full-spectrum',k_f=k_f,k_r=k_r,u_rms=u_rms);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Once a WVTransform and a NonlinearFlux operator have been
% initialized, we can now initialize a model.
% 
% By default the model will use an adaptive time step, although we could
% change this to a fixed time step.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model = WVModel(wvt,nonlinearFlux=fdFlux);

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

model.createNetCDFFileForModelOutput(sprintf('ForcedDissipativeQG-spinup-%d.nc',Nxy),outputInterval=5*86400,shouldOverwriteExisting=1);
% model.setNetCDFOutputVariables('A0','psi','zeta_z','F_psi','F0_psi');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Run the model, plot as the fluid goes unstable
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
tl = tiledlayout(3,3,TileSpacing="compact");
title(tl,'surface vorticity')

t0 = wvt.t;
for i=0:8
    model.integrateToTime(t0 + i*5*86400);
       
    nexttile(tl)
    pcolor(wvt.x/1e3, wvt.y/1e3, wvt.zeta_z(:,:,end).'), shading interp
    colormap("gray")
    title(sprintf('%d days',round(wvt.t/86400)))
    xlabel('km'), ylabel('km')
    xtick([]), ytick([])
    pause(0.1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% At the end of the integration, let's make a plot showing the relative
% vorticity and resulting energy spectrum.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EkT = wvt.transformToRadialWavenumber((wvt.A0_TE_factor/wvt.h) .* (wvt.A0.*conj(wvt.A0)));

figure(name="Forced-dissipative quasigeostrophic turbulence spinup",Position=[100 100 1000 400])
tiledlayout(1,2,TileSpacing="compact")

nexttile
pcolor(wvt.X/1000,wvt.Y/1000,wvt.zeta_z), shading interp
colormap("gray")
xlabel('km'), ylabel('km')
title('relative vorticity')

nexttile
plot(wvt.kRadial,EkT/(wvt.kRadial(2)-wvt.kRadial(1))), xlog, ylog, hold on
plot(wvt.kRadial,model_spectrum2D(wvt.kRadial))
ylabel('m^2/s^2')
xlabel('1/m')
title('total energy spectrum')
vlines([k_f,k_r],'g--')