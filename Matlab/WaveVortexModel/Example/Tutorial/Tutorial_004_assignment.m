%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lx = 2000e3;
Ly = 1000e3;

Nx = 128;
Ny = 64;
Nz = 40;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lz = 4000;
N0 = 3*2*pi/3600; % buoyancy frequency at the surface, radians/seconds
L_gm = 1300; % thermocline exponential scale, meters
N2 = @(z) N0*N0*exp(2*z/L_gm);
wvt = WVTransformHydrostatic([Lx, Ly, Lz], [Nx, Ny, Nz], N2=N2,latitude=25);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Add an eddy
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% "Shallow eddy"
% Density anomaly sits close to the surface
Le = 80e3;
He = wvt.Lz/5;
U = 0.30; % m/s
x0 = (3/4)*max(wvt.x), y0=max(wvt.y)/2;
psi = @(x,y,z) U*(Le/sqrt(2))*exp(1/2)*exp(-((x-x0)/Le).^2 -((y-y0)/Le).^2 -(z/He).^2 );

wvt.setGeostrophicStreamfunction(psi);
% fprintf('min-rv: %.2f f, max-rv: %.2f f\n',min(rv(:))/wvt.f,max(rv(:))/wvt.f)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot SSH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
tl = tiledlayout(3,1);
nexttile

ssh = wvt.seaSurfaceHeight;
pcolor(wvt.x/1e3, wvt.y/1e3, ssh.'), shading interp
axis equal
% max(ssh(:))

rho = wvt.rho_prime;
rho_total = wvt.rho_total;
sliceIndex = find(wvt.y<200e3,1,'last');
sliceIndex = floor(wvt.Ny/2);

nexttile
pcolor(wvt.x/1000,wvt.z,squeeze(rho(:,sliceIndex,:)).'); colorbar; clim([min(rho(:)),max(rho(:))]), shading interp, hold on
contour(wvt.x/1000,wvt.z,squeeze(rho_total(:,sliceIndex,:)).',linspace(min(rho_total(:)),max(rho_total(:)),10),'k','LineWidth',0.5);
xlabel('x (km)'), ylabel('z (m)')

nexttile
rv = wvt.diffX(wvt.v) - wvt.diffY(wvt.u);
pcolor(wvt.x/1000,wvt.z,squeeze(rv(:,sliceIndex,:)).'); colorbar; clim([-1,1]*max(abs(rv(:)))), shading interp
xlabel('x (km)'), ylabel('z (m)')

%%

A0 = wvt.A0;
TE = wvt.A0_TE_factor .* (A0 .* conj(A0));
TE_radial = wvt.transformToRadialWavenumber(TE);
figure, pcolor(wvt.kRadial,wvt.j,TE_radial.'), shading flat

%%
[Fp,Fm,F0] = wvt.nonlinearFlux();

% figure, pcolor(wvt.k,wvt.l,fftshift(TE(:,:,2)).'), shading flat

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up the integrator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize the integrator with the model
model = WVModel(wvt,nonlinearFlux=QGPVE(wvt,shouldUseBeta=1,u_damp=wvt.uMax));
model.setupIntegrator(timeStepConstraint="advective");
% model.createNetCDFFileForModelOutput('qg-eddy.nc')
% model.setNetCDFOutputVariables('u','v','eta','seaSurfaceHeight');
model.integrateToTime(3*86400);

% wvt.writeToFile('eddy-365.nc');
% 
% 
% %%
% 
% wvt = WVTransform.waveVortexTransformFromFile('eddy-365.nc');

U_io = 0.2;
Ld = wvt.Lz/5;
u_NIO = @(z) U_io*exp(-(z/Ld));
v_NIO = @(z) zeros(size(z));

wvt.setInertialMotions(u_NIO,v_NIO);

%%

model = WVModel(wvt,nonlinearFlux=Boussinesq(wvt,shouldAntialias=1,uv_damp=wvt.uMax));
model.setupIntegrator(timeStepConstraint="advective");
model.integrateToTime(wvt.t + 1*wvt.inertialPeriod);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% energy
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
tiledlayout('flow')
Ekj = wvt.transformToRadialWavenumber( wvt.A0_TE_factor .* abs(wvt.A0).^2 );
nexttile, pcolor(wvt.j,wvt.kRadial,Ekj), shading flat
nexttile, plot(wvt.j,sum(Ekj,1))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% y-z slice of (u,v,eta) through the entire domain
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[u,v,eta] = wvt.variables('u','v','eta');
rv = wvt.diffX(wvt.v) - wvt.diffY(wvt.u);

% Let's plot a slice through domain
xIndices = 1:wvt.Nx;
yIndices = floor(wvt.Ny/2);
zIndices = 1:wvt.Nz;
horzAxis = wvt.x/1e3;
vertAxis = wvt.z;

figure
tl = tiledlayout(3,1);
title(tl,sprintf('NIO on day %d',round(wvt.t/86400)))
xlabel(tl,'y (km)')
ylabel(tl,'depth (m)')

nexttile
p1 = pcolor(horzAxis,vertAxis,squeeze(100*u(xIndices,yIndices,zIndices)).'); shading interp
title('u (x-velocity)')
cb1 = colorbar('eastoutside');
cb1.Label.String = 'cm/s';
clim([-5 5])

nexttile
p2 = pcolor(horzAxis,vertAxis,squeeze(rv(xIndices,yIndices,zIndices)).'/wvt.f); shading interp
title('rv')
cb2 = colorbar('eastoutside');
cb2.Label.String = 'f_0';
clim([-0.3 0.3])

nexttile
p3 = pcolor(horzAxis,vertAxis,squeeze(100*eta(xIndices,yIndices,zIndices)).'); shading interp
title('\eta (isopycnal deviation)')
cb3 = colorbar('eastoutside');
cb3.Label.String = 'cm';
clim([-150 150])