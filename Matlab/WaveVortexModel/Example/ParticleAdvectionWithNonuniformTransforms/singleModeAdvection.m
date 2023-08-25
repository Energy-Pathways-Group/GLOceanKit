% set constants to intialize wvt
Lx = 2000e3;
Ly = 1000e3;
Nx = 64;
Ny = 32;
latitude = 25;

x0 = 3*Lx/4;
y0 = Ly/2;
A = 1000;
L = 80e3;
streamFunction = @(x,y) wvt.f/wvt.g * A*exp( - ((x-x0).^2 + (y-y0).^2)/L^2);

% initialize transform and model
wvt = WVTransformSingleMode([Lx, Ly], [Nx, Ny], h=0.8, latitude=latitude);
wvt.setSSH(streamFunction);
model = WVModel(wvt,nonlinearFlux=SingleModeQGPVE(wvt,shouldUseBeta=1,u_damp=0));

% set drifter positions
xFloatRange = linspace(x0, x0+Lx/30, Nx/4);
yFloatRange = linspace(y0, y0+Ly/30, Ny/4);
[xFloat,yFloat] = ndgrid(xFloatRange,yFloatRange);
xFloat = reshape(xFloat,1,[]);
yFloat = reshape(yFloat,1,[]);
nTrajectories = length(xFloat);
model.setDrifterPositions(xFloat,yFloat,[],"ssh",...
    advectionInterpolation="finufft",...
    trackedVarInterpolation="finufft"...
    );


% graph fluid velocity and particle locations
colorbar_lims = [-0.01, 0.01];
figure(2);
clf;
hold on;
pcolor(...
        wvt.x, ...
        wvt.y, ...
        wvt.ssh.' ...
        );
scatter(xFloat, yFloat, '.');
colorbar;
clim(colorbar_lims);
shading interp;
title("sea surface height (m) at time 0");

figure(3);
clf;
hold on;
pcolor(...
        wvt.x, ...
        wvt.y, ...
        wvt.u.' ...
        );
scatter(xFloat, yFloat, '.');
colorbar;
clim(colorbar_lims);
shading interp;
title("u (m/s) at time 0");

figure(4);
clf;
hold on;
pcolor(...
        wvt.x, ...
        wvt.y, ...
        wvt.v.' ...
        );
scatter(xFloat, yFloat, '.');
colorbar;
clim(colorbar_lims);
shading interp;
title("v (m/s) at time 0");


% prepare output
period=wvt.inertialPeriod;
model.setupIntegrator(timeStepConstraint="advective", outputInterval=period/10);
model.createNetCDFFileForModelOutput('AccurateParticleAdvection.nc', ...
    shouldOverwriteExisting=1);

% integrate
finalperiodcount = 40;
model.integrateToTime(finalperiodcount*period);


figure(5);
pcolor(...
        wvt.x, ...
        wvt.y, ...
        wvt.ssh.' ...
        );
shading interp;
colorbar;
clim(colorbar_lims);
title("sea surface height (m) after "+finalperiodcount+" periods");

figure(6);
pcolor(...
        wvt.x, ...
        wvt.y, ...
        wvt.u.' ...
        );
shading interp;
colorbar;
clim(colorbar_lims);
title("u (m/s) after "+finalperiodcount+" periods");

figure(7);
pcolor(...
        wvt.x, ...
        wvt.y, ...
        wvt.v.' ...
        );
shading interp;
colorbar;
clim(colorbar_lims);
title("v (m/s) after "+finalperiodcount+" periods");


%% load data from file
ncfile = model.ncfile;
[x,y] = ncfile.readVariables('drifter-x','drifter-y');

figure(1);
clf;
hold on;
scatter(x(:,1), y(:,1));
plot(x.',y.');

