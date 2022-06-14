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

wvt = WaveVortexTransformSingleMode([Lx, Ly], [Nx, Ny], h=0.8, latitude=latitude);

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
model = WaveVortexModel(wvt,nonlinearFlux=SingleModeQGPVE(wvt,shouldUseBeta=1));
% model.nonlinearFlux = SingleModeQGPVE(model.wvt,shouldUseBeta=1);


% set initial positions for a bunch of floats
nTrajectories = 15000;
x = linspace(Lx/8,(7/8)*Lx, round((Lx/Ly)*sqrt(nTrajectories/(Lx/Ly))) );
y = linspace(Ly/8,(7/8)*Ly, round(sqrt(nTrajectories/(Lx/Ly))));
[xFloat,yFloat] = ndgrid(x,y);
xFloat = reshape(xFloat,1,[]);
yFloat = reshape(yFloat,1,[]);
nTrajectories = length(xFloat);
model.setDrifterPositions(xFloat,yFloat,[],'qgpv');

[deltaT,advectiveDT,oscillatoryDT] = model.timeStepForCFL(0.15);
finalTime = 75*86400;
nT = model.setupIntegrator(advectiveDT, 86400, finalTime);

model.createNetCDFFileForModelOutput('QGMonopole.nc',shouldOverwriteExisting=1);
model.integrateToTime(finalTime);

ncfile = model.ncfile;
[x,y] = ncfile.readVariables('drifter-x','drifter-y');

figure, plot(x.',y.')
