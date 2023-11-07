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

% "Deep eddy"
% We are going to have its velocity vanish at the bottom, and only have its
% density anomaly in the middle--hence the errorfunction.
x0 = 3*Lx/4;
y0 = Ly/2;
% Le = 80e3;
% z0 = -wvt.Lz/4;
% He = wvt.Lz/10;
% U = 0.25; % m/s
% psi = @(x,y,z) U*(Le/sqrt(2))*exp(1/2)*exp(-((x-x0)/Le).^2 -((y-y0)/Le).^2) .* (erf((z-z0)/He)+1)/2;

% "Shallow eddy"
% Density anomaly sits close to the surface
Le = 80e3;
He = wvt.Lz/5;
U = 0.30; % m/s
psi = @(x,y,z) U*(Le/sqrt(2))*exp(1/2)*exp(-((x-x0)/Le).^2 -((y-y0)/Le).^2 -(z/He).^2 );

wvt.setGeostrophicStreamfunction(psi);
% fprintf('min-rv: %.2f f, max-rv: %.2f f\n',min(rv(:))/wvt.f,max(rv(:))/wvt.f)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot SSH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ssh = wvt.seaSurfaceHeight;
figure, pcolor(wvt.x/1e3, wvt.y/1e3, ssh.'), shading interp
max(ssh(:))

% figure, pcolor(wvt.x,wvt.y,wvt.ssh.'), shading interp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize the integrator with the model
model = WVModel(wvt,nonlinearFlux=WVNonlinearFluxQG(wvt,shouldUseBeta=1,uv_damp=wvt.uMax,r = 1/(200*86400)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add floats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% xf = wvt.x(1:4:end);
% yf = wvt.y(1:4:end);
% zf = wvt.z(end:-4:1); % make sure we grab the surface
% [XF,YF,ZF] = ndgrid(xf,yf,zf);
% model.setFloatPositions(XF(:),YF(:),ZF(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up the integrator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model.setupIntegrator(timeStepConstraint="advective",outputInterval=86400);
model.createNetCDFFileForModelOutput('qg-eddy.nc',shouldOverwriteExisting=1,shouldUseClassicNetCDF=1);
model.setNetCDFOutputVariables('u','v','eta','rho_prime','seaSurfaceHeight');
model.integrateToTime(100*86400);

ssh = wvt.seaSurfaceHeight;
figure, pcolor(wvt.x/1e3, wvt.y/1e3, ssh.'), shading interp
max(ssh(:))

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read float positions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ncfile = model.ncfile;
% [x,y,z] = ncfile.floatPositions();
[xFloat,yFloat,zFloat] = ncfile.readVariables('float_x','float_y','float_z');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uniqueDepths = unique(zFloat(:,1));

figure
tl = tiledlayout('flow');
tl.TileSpacing = 'compact';
tl.Padding = 'compact';

for iDepth=length(uniqueDepths):-1:1
    indices = find( abs(zFloat(:,1)-uniqueDepths(iDepth)) < 1e-3);
    x = xFloat(indices,:);
    y = yFloat(indices,:);
    z = zFloat(indices,:);
    nexttile
    plot(x.',y.')
    title(sprintf('float at depth %d m',round(uniqueDepths(iDepth))))
end

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