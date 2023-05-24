%% Initialize a wave vortex transform

% setup stratification
N0 = 3*2*pi/3600;
L_gm = 1300;
N2 = @(z) N0*N0*exp(2*z/L_gm);
wvt = WVTransformHydrostatic([2000e3, 1000e3, 4000],[64, 32, 32], N2=N2, latitude=30);

% setup initial flow
Le = 80e3;
He = wvt.Lz/5;
U = 0.263; % m/s
x0 = (3/4)*max(wvt.x)
y0 = max(wvt.y)/2
psi = @(x,y,z) U*(Le/sqrt(2))*exp(1/2)*exp(-((x-x0)/Le).^2 -((y-y0)/Le).^2 -(z/He).^2 );
wvt.initWithGeostrophicStreamfunction(psi);


%% Estimate energy content in modes

% wvt.summarizeEnergyContent

wvt.summarizeModeEnergy

% wvt.totalEnergy
% 
% wvt.totalEnergySpatiallyIntegrated
% 
% wvt.totalHydrostaticEnergy
% 
% wvt.geostrophicEnergy
% 
% wvt.inertialEnergy
% 
% wvt.waveEnergy


%% Make some plots

figure(1);clf
pcolor(wvt.x/1000, wvt.y/1000, wvt.seaSurfaceHeight'); shading flat
xlabel('x (km)')
ylabel('y (km)')
c=colorbar;
c.Label.String = 'SSH (m)';
caxis([0,.15])

figure(2);clf
rho_total = wvt.rho_total;
pcolor(wvt.x/1000, wvt.z, squeeze(rho_total(:,wvt.Ny/2,:))'); shading flat
xlabel('x (km)')
ylabel('z (m)')
c=colorbar;
c.Label.String = 'rho total (kg m^{-3})';

figure(3);clf
rho_prime = wvt.rho_prime;
pcolor(wvt.x/1000, wvt.z, squeeze(rho_prime(:,wvt.Ny/2,:))'); shading flat
xlabel('x (km)')
ylabel('z (m)')
c=colorbar;
c.Label.String = 'rho total (kg m^{-3})';
cmocean('Balance');
caxis([-max(abs(caxis)), max(abs(caxis))])


%% Compute the vertical component of the relative vorticity (partial_x v - partial_y u) and plot

rv = wvt.diffX(wvt.v) - wvt.diffY(wvt.u);

figure(4);clf
pcolor(wvt.x/1000, wvt.y/1000, squeeze(rv(:,:,floor(wvt.Nz/3)))'); shading flat
xlabel('x (km)')
ylabel('y (km)')
c=colorbar;
c.Label.String = 'relative vorticity (s^{-1})';
cmocean('Balance');


%% Initialize a model and use the QGPVE --- on the beta plane!

model = WVModel(wvt, nonlinearFlux=QGPVE(wvt,shouldUseBeta=1,u_damp=wvt.uMax));


%% Add particles with float like behavior


%% Run for 365 days

% setup integrator
model.setupIntegrator(timeStepConstraint="advective"); % outputInterval=wvt.inertialPeriod/2
% create output file
model.createNetCDFFileForModelOutput('Tutorial_day3practice.nc',shouldOverwriteExisting=1);
% integrate for 365 days
model.integrateToTime(365*86400);


%% Plot movement of particle



