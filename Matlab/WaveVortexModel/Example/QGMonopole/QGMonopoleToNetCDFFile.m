%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lx = 2000e3;
Ly = 1000e3;

Nx = 256;
Ny = 128;

latitude = 25;

wvt = WVTransformSingleMode([Lx, Ly], [Nx, Ny], h=0.8, latitude=latitude);

x0 = 3*Lx/4;
y0 = Ly/2;
A = 0.15;
L = 80e3;
wvt.setSSH(@(x,y) A*exp( - ((x-x0).^2 + (y-y0).^2)/L^2) );

figure, pcolor(wvt.x,wvt.y,wvt.ssh.'), shading interp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the integrator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize the integrator with the model
model = WVModel(wvt,nonlinearFlux=QGPVE(wvt,shouldUseBeta=1,u_damp=wvt.uMax));

% set initial positions for a bunch of floats
[xFloat,yFloat] = ndgrid(wvt.x(1:2:end),wvt.y(1:2:end));
xFloat = reshape(xFloat,1,[]);
yFloat = reshape(yFloat,1,[]);
nTrajectories = length(xFloat);
model.setDrifterPositions(xFloat,yFloat,[],'qgpv');

model.setupIntegrator(timeStepConstraint="advective", outputInterval=86400);

model.createNetCDFFileForModelOutput('QGMonopole.nc',shouldOverwriteExisting=1);
model.setNetCDFOutputVariables('A0');
model.integrateToTime(75*86400);

ncfile = model.ncfile;
[x,y] = ncfile.readVariables('drifter_x','drifter_y');
qgpv = ncfile.readVariables('drifter_qgpv');

figure, plot(x.',y.')
