Lx = 750e3;
Ly = 750e3;

N = 128;
Nx = N;
Ny = N;
nModes = 40;

latitude = 31;

% Simulation length
inertialPeriod = (2*pi/(2 * 7.2921E-5 * sin( latitude*pi/180 )));
maxTime = 1*inertialPeriod;
outputInterval = inertialPeriod/10;

outputfile = '/Volumes/MoreStorage/Data/cyprus_eddy_wvm/cyprus_eddy-9.nc';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Setup the model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N0 = 12*2*pi/3600; % reference buoyancy frequency, radians/seconds
Nmin = N0*3e-2;
rho0 = 1025; g = 9.81;
L_gm = 145; % thermocline exponential scale, meters
rhoFunction = @(z) rho0*(1 + L_gm*N0*N0/(2*g)*(1 - exp(2*z/L_gm)) - (Nmin*Nmin/g)*z);
N2Function = @(z) N0*N0*exp(2*z/L_gm) + Nmin*Nmin*ones(size(z));
dLnN2Function = @(z) 2*ones(size(z))/L_gm;

% save('/Volumes/MoreStorage/Data/cyprus_eddy_wvm/cyprus_eddy-2.mat','rhoFunction','N2Function','dLnN2Function');

zIn = [-2000 0];

wvm = WaveVortexModelHydrostatic([Lx, Ly, max(zIn)-min(zIn)], [Nx, Ny, nModes], latitude, rhoFunction,'N2func', N2Function, 'dLnN2func',dLnN2Function);

% figure
% 
% minN2 = 0;
% maxN2 = 1.1*max(sqrt(wvm.N2))*3600/2/pi; % radians/sec * (1/2pi)*(86400)
% xh = repmat([minN2 maxN2],length(wvm.z),1);
% zh = [wvm.z,wvm.z];
% 
% sp1 = subplot(2,1,1);
% plot(sqrt(abs(wvm.N2))*3600/2/pi,wvm.z,'LineWidth',2,'Color',0.2*[1 1 1]), hold on
% plot(xh.',zh.','k')
% scatter(interp1(wvm.z,sqrt(abs(wvm.N2))*3600/2/pi,wvm.z),wvm.z,5^2,0.0*[1 1 1],'filled')
% xlim([minN2 maxN2])
% ylim([-150 0])
% yticks(linspace(-100,0,3))
% % set( gca, 'FontSize', figure_axis_tick_size);
% 
% sp2 = subplot(2,1,2);
% plot(sqrt(abs(wvm.N2))*3600/2/pi,wvm.z,'LineWidth',2,'Color',0.2*[1 1 1]), hold on
% plot(xh.',zh.','k')
% scatter(interp1(wvm.z,sqrt(abs(wvm.N2))*3600/2/pi,wvm.z),wvm.z,5^2,0.0*[1 1 1],'filled')
% xlim([minN2 maxN2])
% ylim([-1650 -150])
% yticks(linspace(-1650,-150,3))
% xlabel('buoyancy frequency (cph)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
% ylabel('depth (m)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
% set( gca, 'FontSize', figure_axis_tick_size);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Add a geostrophic stream function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x0 = (max(wvm.x)-min(wvm.x))/2;
y0 = (max(wvm.y)-min(wvm.y))/2;

A = -1.32e-5; % s^{-1}
alpha = 8e-10; % m^{-2}
beta = 8.2e-6; % m^{-2}
psi = @(x,y,z) -(A/(2*alpha))*exp(-alpha*((x-x0).*(x-x0)+(y-y0).*(y-y0))-beta*z.*z);

u = @(x,y,z) -A*(y-y0).*exp(-alpha*((x-x0).*(x-x0)+(y-y0).*(y-y0))-beta*z.*z);
v = @(x,y,z) A*(x-x0).*exp(-alpha*((x-x0).*(x-x0)+(y-y0).*(y-y0))-beta*z.*z);
% rho = @(x,y,z) -(rho0/g)*A*wvm.f0*(beta/alpha)*z.*exp(-alpha*((x-x0).*(x-x0)+(y-y0).*(y-y0))-beta*z.*z) - (rho0/g)*A*A*(beta/alpha)*z.*exp(-2*alpha*((x-x0).*(x-x0)+(y-y0).*(y-y0))-2*beta*z.*z);
eta = @(x,y,z) -A*wvm.f0*(beta/alpha)*(z./N2Function(z)).*exp(-alpha*((x-x0).*(x-x0)+(y-y0).*(y-y0))-beta*z.*z) - A*A*(beta/alpha)*(z./N2Function(z)).*exp(-2*alpha*((x-x0).*(x-x0)+(y-y0).*(y-y0))-2*beta*z.*z);
[X,Y,Z]=ndgrid(wvm.x,wvm.y,wvm.z);
U = u(X,Y,Z);
V = v(X,Y,Z);
N = eta(X,Y,Z);

wvm.shouldAntiAlias = 1;
[Qk,Ql,Qj] = wvm.ExponentialFilter();
Q = Qk.*Ql.*Qj;
[Ap,Am,A0] = wvm.TransformUVEtaToWaveVortex(U,V,N);
wvm.Ap = Q.*Ap;
wvm.Am = Q.*Am;
wvm.A0 = Q.*A0;

% [wvm.Ap,wvm.Am,wvm.A0] = wvm.TransformUVEtaToWaveVortex(U,V,N);


rho = wvm.DensityFieldAtTime(0);
figure, plot(squeeze(rho(Nx/2,Ny/2,:)),wvm.z)
hold on, plot(wvm.rhobar,wvm.z)

[rho_prime,eta] = wvm.VariableFieldsAtTime(0,'rho_prime','eta');
figure, plot(squeeze(N(Nx/2,Ny/2,:)),wvm.z), hold on
plot(squeeze(eta(Nx/2,Ny/2,:)),wvm.z)
figure, plot(squeeze(rho_prime(Nx/2,Ny/2,:)),wvm.z), hold on

dStrat = log10(max(wvm.N2)/min(wvm.N2));
if (dStrat > 7)
    warning('Mean stratification (N2) changes by %d orders of magnitude. This may lead to numerical instability.',round(dStrat));
end

wvm.totalHydrostaticEnergy
wvm.totalSpectralEnergy

return


% wvm.SetGeostrophicStreamfunction(psi);



% x = wvm.x;
% y = wvm.y;
% z = wvm.z;
% f0 = wvm.f0;
% zeta_z = (DiffFourier(x,v,1,1) - DiffFourier(y,u,1,2))/f0;
% figure, pcolor(x/1e3,z/1e3,squeeze(zeta_z(:,Ny/2,:)).'), shading interp; colorbar('eastoutside')
% figure, pcolor(x/1e3,y/1e3,squeeze(zeta_z(:,:,35)).'), shading interp; colorbar('eastoutside')
% % figure, pcolor(x/1e3,y/1e3,squeeze(zeta_z(:,:,35)).'), shading interp; colorbar('eastoutside')
% 
% var = eta-N;
% figure, pcolor(x/1e3,z/1e3,squeeze(var(:,Ny/2,:)).'), shading interp; colorbar('eastoutside')
% figure, pcolor(x/1e3,z/1e3,squeeze(N(:,Ny/2,:)).'), shading interp; colorbar('eastoutside')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Add a an inertial function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% B = 0.2; % m/s
% gamma = 0.01; % m^{-1}
% u_IO = @(z) B*exp(gamma*z);
% v_IO = @(z) zeros(size(z));
% 
% wvm.SetInertialMotions(u_IO,v_IO);

% N=5;tmp=cat(1,squeeze(Ubar(1,1,1:end-N)),zeros(N,1));figure, plot(self.PFinv * tmp,self.z)

% [u,v,eta] = wvm.VariableFieldsAtTime(0, 'u','v','eta');
% figure, pcolor(x/1e3,z/1e3,squeeze(u(:,Ny/2,:)).'), shading interp; colorbar('eastoutside')
% figure, pcolor(x/1e3,y/1e3,squeeze(u(:,:,floor(0.9*nModes))).'), shading interp; colorbar('eastoutside')

wvm.summarizeEnergyContent;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Determine the proper time interval
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

omega = wvm.Omega;
period = 2*pi/max(abs(omega(:)));
[u,v] = wvm.VelocityFieldAtTime(0.0);
U = max(max(max( sqrt(u.*u + v.*v) )));
fprintf('Max fluid velocity: %.2f cm/s\n',U*100);

cfl = 0.25;
advectiveDT = cfl*(wvm.x(2)-wvm.x(1))/U;
oscillatoryDT = period/12;
if advectiveDT < oscillatoryDT
    fprintf('Using the advective dt: %.2f s\n',advectiveDT);
    deltaT = advectiveDT;
else
    fprintf('The highest frequency resolved IGW has period %.1f min.\n',period/60);
    fprintf('Using the oscillatory dt: %.2f s\n',oscillatoryDT);
    deltaT = oscillatoryDT;
end

deltaT = outputInterval/ceil(outputInterval/deltaT);
stepsPerOutput = round(outputInterval/deltaT);
fprintf('Rounding to match the output interval dt: %.2f s (%d steps per output)\n',deltaT,stepsPerOutput);

t = (0:outputInterval:maxTime)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Setup a netcdf file for the output
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

netcdfTool = WaveVortexModelNetCDFTools(outputfile);
netcdfTool.CreateNetCDFFileFromModel(wvm,length(t),'double');
netcdfTool.CreateAmplitudeCoefficientVariables();
netcdfTool.CreateEnergeticsVariables();
netcdfTool.CreateEnergeticsKJVariables();

% Save the initial conditions
iTime = 1;
netcdfTool.WriteTimeAtIndex(iTime,t(iTime));
netcdfTool.WriteAmplitudeCoefficientsAtIndex(iTime);
netcdfTool.WriteEnergeticsAtIndex(iTime);
netcdfTool.WriteEnergeticsKJAtIndex(iTime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Time step forward!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

integrator = ArrayIntegrator(@(t,y0) wvm.NonlinearFluxAtTimeArray(t,y0),{wvm.Ap,wvm.Am,wvm.A0},deltaT);

% profile on
startTime = datetime('now');
fprintf('Starting numerical simulation on %s\n', datestr(startTime));
for iTime=2:length(t)
    if iTime == 2
        startTime = datetime('now');
    end
    if iTime == 3 || mod(iTime,5) == 0
        timePerStep = (datetime('now')-startTime)/(iTime-2);
        timeRemaining = (length(t)-iTime+1)*timePerStep;
        fprintf('\twriting values time step %d of %d to file. Estimated finish time %s (%s from now)\n', iTime, length(t), datestr(datetime('now')+timeRemaining), datestr(timeRemaining, 'HH:MM:SS')) ;
    end

    for i=1:stepsPerOutput
        integrator.IncrementForward();
        wvm.Ap = integrator.currentY{1};
        wvm.Am = integrator.currentY{2};
        wvm.A0 = integrator.currentY{3};
    end

    if mod(iTime,10)==0
        wvm.summarizeEnergyContent();
    end
    netcdfTool.WriteTimeAtIndex(iTime,t(iTime));
    netcdfTool.WriteAmplitudeCoefficientsAtIndex(iTime);
    netcdfTool.WriteEnergeticsAtIndex(iTime);
    netcdfTool.WriteEnergeticsKJAtIndex(iTime);

[Ep,Em,E0] = wvm.EnergyFluxAtTime(integrator.currentTime,wvm.Ap,wvm.Am,wvm.A0);
inertialFlux = sum(Ep(1,1,:)) + sum(Em(1,1,:));
Ep(1,1,:) = 0; Em(1,1,:) = 0;
waveFlux = sum(Ep(:)) + sum(Em(:));
fprintf('total spectral: %g, total depth integrated: %g\n',wvm.totalSpectralEnergy,wvm.totalHydrostaticEnergy);
fprintf('io flux: %g, wave flux: %g, geostrophic flux: %g. Net: %g\n',inertialFlux,waveFlux,sum(E0(:)),inertialFlux+waveFlux+sum(E0(:)))

% if iTime>170
%     [u,v] = wvm.VelocityFieldAtTime(iTime);
%     U = max(max(max( sqrt(u.*u + v.*v) )));
%     fprintf('cfl%.2f\n',(deltaT*U/wvm.x(2)-wvm.x(1)));
%     wvm.summarizeEnergyContent();
% end
end

[k,j,ApKJ,AmKJ,A0KJ] = wvm.ConvertToWavenumberAndMode(abs(wvm.Ap).^2,abs(wvm.Am).^2,abs(wvm.A0).^2);
% [k,j,ApKJ,AmKJ,A0KJ] = self.ConvertToWavenumberAndMode(abs(uNLbar.^2),abs(vNLbar).^2,abs(nNLbar).^2);
[k,j,ApKJ,AmKJ,A0KJ] = wvm.ConvertToWavenumberAndMode(Ep,Em,E0);
figure
subplot(2,2,1)
plot(k,sum(ApKJ,2)), hold on, plot(k,sum(AmKJ,2))%, ylog
subplot(2,2,2)
plot(k,sum(A0KJ,2))%, ylog

subplot(2,2,3)
plot(j,sum(ApKJ,1)), hold on, plot(j,sum(AmKJ,1))%, ylog
subplot(2,2,4)
plot(j,sum(A0KJ,1))%, ylog
