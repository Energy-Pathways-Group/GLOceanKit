%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wvt = WVTransformConstantStratification([50e3, 50e3, 1300],[64, 64, 65], N0=5.2e-3,latitude=25);
omega = wvt.initWithWaveModes(k=10,l=0,j=1,phi=0,U=0.2,sign=1);
period = 2*pi/omega;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the integrator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize the integrator with the model
model = WVModel(wvt,nonlinearFlux=WVNonlinearFlux(wvt,shouldAntialias=1,uv_damp=wvt.uMax));

% set initial positions for a bunch of floats
nTrajectories = 101;
xFloat = wvt.Lx/2*ones(1,nTrajectories);
yFloat = wvt.Ly/2*ones(1,nTrajectories);
zFloat = linspace(-wvt.Lz,0,nTrajectories);
model.setFloatPositions(xFloat,yFloat,zFloat);

model.setupIntegrator(timeStepConstraint="oscillatory", outputInterval=period/10);
model.createNetCDFFileForModelOutput('PlaneWaveWithFloats.nc',shouldOverwriteExisting=1);
model.integrateToTime(3*period);

ncfile = model.ncfile;
% [x,y,z] = ncfile.floatPositions();
[x,z] = ncfile.readVariables('float_x','float_z');

figure, plot(x.',z.')