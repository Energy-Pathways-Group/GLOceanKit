Lx = 750e3;
Ly = 750e3;

N = 128;
Nx = N;
Ny = N;
nModes = 60;

latitude = 33.5;

% Simulation length
inertialPeriod = (2*pi/(2 * 7.2921E-5 * sin( latitude*pi/180 )));
maxTime = 2*inertialPeriod;
outputInterval = inertialPeriod/10;

outputfile = '/Volumes/MoreStorage/Data/cyprus_eddy_wvm/cyprus_eddy-more-stratification-strong-2.nc';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Setup the model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N0 = 12*2*pi/3600; % buoyancy frequency at the surface, radians/seconds
rho0 = 1025; g = 9.81;
L_gm = 2*145; % thermocline exponential scale, meters
L_const = 3*L_gm; % depth below which stratification stays constant—***Note*** if you use 8*L_gm, stratification will be too weak compared to the size of the buoyancy anomaly for sensible results.
Nmin = sqrt(N0*N0*exp(-2*L_const/L_gm));
rhoFunction = @(z) rho0*(1 + L_gm*N0*N0/(2*g)*(1 - exp(2*z/L_gm)) - (Nmin*Nmin/g)*z);
N2Function = @(z) N0*N0*exp(2*z/L_gm) + Nmin*Nmin*ones(size(z));
dLnN2Function = @(z) 2*ones(size(z))/L_gm;

zIn = [-2000 0];

wvm = WaveVortexModelHydrostatic([Lx, Ly, max(zIn)-min(zIn)], [Nx, Ny, nModes], latitude, rhoFunction,'N2func', N2Function, 'dLnN2func',dLnN2Function);
wvm.shouldAntiAlias = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the cyclo-geostrophic flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x0 = (max(wvm.x)-min(wvm.x))/2;
y0 = (max(wvm.y)-min(wvm.y))/2;

A = -1.32e-5; % s^{-1}
alpha = 1/((35e3)^2); % m^{-2}
beta = 1/(350^2); % m^{-2}
psi = @(x,y,z) -(A/(2*alpha))*exp(-alpha*((x-x0).*(x-x0)+(y-y0).*(y-y0))-beta*z.*z);

u_eddy = @(x,y,z) -A*(y-y0).*exp(-alpha*((x-x0).*(x-x0)+(y-y0).*(y-y0))-beta*z.*z);
v_eddy = @(x,y,z) A*(x-x0).*exp(-alpha*((x-x0).*(x-x0)+(y-y0).*(y-y0))-beta*z.*z);
% rho_eddy = @(x,y,z) -(rho0/g)*A*wvm.f0*(beta/alpha)*z.*exp(-alpha*((x-x0).*(x-x0)+(y-y0).*(y-y0))-beta*z.*z) - (rho0/g)*A*A*(beta/alpha)*z.*exp(-2*alpha*((x-x0).*(x-x0)+(y-y0).*(y-y0))-2*beta*z.*z);
eta_eddy = @(x,y,z) -A*wvm.f0*(beta/alpha)*(z./N2Function(z)).*exp(-alpha*((x-x0).*(x-x0)+(y-y0).*(y-y0))-beta*z.*z) - A*A*(beta/alpha)*(z./N2Function(z)).*exp(-2*alpha*((x-x0).*(x-x0)+(y-y0).*(y-y0))-2*beta*z.*z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the inertial flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B = 0.4; % m/s
gamma = 1/100; % m^{-1}
u_IO = @(z) B*exp(gamma*z);
v_IO = @(z) zeros(size(z));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add the two components together and set as initial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X,Y,Z]=ndgrid(wvm.x,wvm.y,wvm.z);

U = u_eddy(X,Y,Z) + u_IO(Z);
V = v_eddy(X,Y,Z) + v_IO(Z);
N = eta_eddy(X,Y,Z);

% Because these are smooth initial conditions, we do need to apply a
% spectral filter to prevent large oscillations in the density anomaly.
nFilteredModes = 20;
[Qk,Ql,Qj] = wvm.ExponentialFilter(nFilteredModes);
Q = Qk.*Ql.*Qj;
[Ap,Am,A0] = wvm.transformUVEtaToWaveVortex(U,V,N);
wvm.Ap = Q.*Ap;
wvm.Am = Q.*Am;
wvm.A0 = Q.*A0;

% Alternative methods for setting initial conditions that I may use in the
% future.
% wvm.setGeostrophicStreamfunction(psi);
% wvm.setInertialMotions(u_IO,v_IO);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sanity checks
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxEta = max(abs(N(:)));
if (maxEta > wvm.Lz/2 )
    warning('Max isopycnal displacement is %.0f meters, compared to a total depth of %0.f meters. This will work, but spectral energies will not be reliable.',maxEta, wvm.Lz);
end
dStrat = log10(max(wvm.N2)/min(wvm.N2));
if (dStrat > 7)
    warning('Mean stratification (N2) changes by %d orders of magnitude. This may lead to numerical instability.',round(dStrat));
end
dEnergy = abs(1- wvm.totalHydrostaticEnergy/wvm.totalEnergy);
if (dEnergy > 0.05)
    warning('Depth integrated and spectral energy disagree by %.1 percent, suggesting misspecified mean density.',dEnergy*100);
end

% Create a figure looking a various aspects of density
rho = wvm.DensityFieldAtTime(0);
[u,v,eta,rho_prime] = wvm.VariableFieldsAtTime(0,'u','v','eta', 'rho_prime');

figure
subplot(2,2,[1 2])
plot(squeeze(rho(Nx/2,Ny/2,:)),wvm.z)
hold on, plot(wvm.rhobar,wvm.z)
subplot(2,2,3)
plot(squeeze(N(Nx/2,Ny/2,:)),wvm.z), hold on
plot(squeeze(eta(Nx/2,Ny/2,:)),wvm.z)
xlabel('eta (m)'), ylabel('depth (m)'), legend('projected, filtered', 'original')
subplot(2,2,4)
plot(squeeze(rho_prime(Nx/2,Ny/2,:)),wvm.z)
xlabel('density anomaly')

% Create a figure looking at the vertical vorticity
zeta_z = (DiffFourier(wvm.x,v,1,1) - DiffFourier(wvm.y,u,1,2))/wvm.f0;
figure('Position',[100 100 400 800])

subplot(2,1,1)
pcolor(wvm.x/1e3,wvm.z,squeeze(zeta_z(:,Ny/2,:)).'), shading interp; colorbar('eastoutside')
xlabel('x (km)'), ylabel('z (m)'), title(sprintf('vorticity at y=%.1f km',wvm.y(Ny/2)/1e3))
subplot(2,1,2)
pcolor(wvm.x/1e3,wvm.y/1e3,squeeze(zeta_z(:,:,35)).'), shading interp; colorbar('eastoutside') 
xlabel('x (km)'), ylabel('y (km)'), title(sprintf('vorticity at z=%.1f m',wvm.z(35)))

wvm.summarizeEnergyContent;

pause(1)

nctool = WaveVortexModelNetCDFTools( outputfile );

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize new integrator, overwrite existing file, integrate
%
integrationTool = WaveVortexModelIntegrationTools(wvm, outputfile, outputInterval,'shouldOverwriteExisting',1);
integrationTool.integrateToTime(2*outputInterval);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Append/continue integration to the same file
%
integrationTool.integrateToTime(4*outputInterval);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Restart from existing output, create new file
%
existingModelOutput = outputfile;
restartModelOutput = '/Volumes/MoreStorage/Data/cyprus_eddy_wvm/cyprus_eddy-more-stratification-strong-2-restart.nc';
integrationTool = WaveVortexModelIntegrationTools(existingModelOutput, restartModelOutput);
integrationTool.integrateToTime(6*outputInterval);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Restart from existing output, create new file, double the resolution
%
restartX2ModelOutput = '/Volumes/MoreStorage/Data/cyprus_eddy_wvm/cyprus_eddy-more-stratification-strong-2-restart-x2.nc';
integrationTool = WaveVortexModelIntegrationTools(existingModelOutput, restartX2ModelOutput,[],'shouldDoubleResolution',1);
integrationTool.integrateToTime(6*outputInterval);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Restart from existing output, append to existing output
%
% integrationTool = WaveVortexModelIntegrationTools(newModelOutput, newModelOutput);
% integrationTool.integrateToTime(4*maxTime/8);

netcdfTool = WaveVortexModelNetCDFTools(restartModelOutput,'timeIndex',Inf);
wvm = netcdfTool.wvm;
t = netcdfTool.t;
wvm = wvm.waveVortexModelWithResolution(2*[wvm.Nx,wvm.Ny,wvm.nModes]);
[u,v] = wvm.VariableFieldsAtTime(t,'u','v');
zeta_z = (DiffFourier(wvm.x,v,1,1) - DiffFourier(wvm.y,u,1,2))/wvm.f0;
figure('Position',[100 100 400 800])
subplot(2,1,1)
pcolor(wvm.x/1e3,wvm.z,squeeze(zeta_z(:,wvm.Ny/2,:)).'), shading interp; colorbar('eastoutside')
xlabel('x (km)'), ylabel('z (m)'), title(sprintf('vorticity at y=%.1f km',wvm.y(wvm.Ny/2)/1e3))
subplot(2,1,2)
pcolor(wvm.x/1e3,wvm.y/1e3,squeeze(zeta_z(:,:,35)).'), shading interp; colorbar('eastoutside') 
xlabel('x (km)'), ylabel('y (km)'), title(sprintf('vorticity at z=%.1f m',wvm.z(35)))

netcdfToolX2 = WaveVortexModelNetCDFTools(restartX2ModelOutput,'timeIndex',Inf);
wvm = netcdfToolX2.wvm;
t = netcdfTool.t;
[u,v] = wvm.VariableFieldsAtTime(t,'u','v');
zeta_z = (DiffFourier(wvm.x,v,1,1) - DiffFourier(wvm.y,u,1,2))/wvm.f0;
figure('Position',[100 100 400 800])
subplot(2,1,1)
pcolor(wvm.x/1e3,wvm.z,squeeze(zeta_z(:,wvm.Ny/2,:)).'), shading interp; colorbar('eastoutside')
xlabel('x (km)'), ylabel('z (m)'), title(sprintf('vorticity at y=%.1f km',wvm.y(wvm.Ny/2)/1e3))
subplot(2,1,2)
pcolor(wvm.x/1e3,wvm.y/1e3,squeeze(zeta_z(:,:,35)).'), shading interp; colorbar('eastoutside') 
xlabel('x (km)'), ylabel('y (km)'), title(sprintf('vorticity at z=%.1f m',wvm.z(35)))

