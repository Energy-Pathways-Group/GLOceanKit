%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 64;

Lx = 50e3;
Ly = Lx;
Lz = 1300;

Nx = N;
Ny = N;
Nz = N+1; % 2^n + 1 grid points, to match the Winters model, but 2^n ok too.

latitude = 25;
N0 = 5.2e-3; % Choose your stratification 7.6001e-04

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wvt = WaveVortexTransformConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], N0,latitude=latitude);

U = .2; phi=0;
omega = wvt.initWithWaveModes(10,0,1,phi,U,1);
period = 2*pi/omega;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the integrator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize the integrator with the model
model = WVModel(wvt);

% set initial positions for a bunch of floats
nTrajectories = 101;
xFloat = Lx/2*ones(1,nTrajectories);
yFloat = Ly/2*ones(1,nTrajectories);
zFloat = linspace(-Lz,0,nTrajectories);
model.setFloatPositions(xFloat,yFloat,zFloat);

model.setupIntegrator(timeStepConstraint="oscillatory", outputInterval=period/10);
model.createNetCDFFileForModelOutput('PlaneWaveWithFloats.nc',shouldOverwriteExisting=1);
model.integrateToTime(3*period);

ncfile = model.ncfile;
% [x,y,z] = ncfile.floatPositions();
[x,z] = ncfile.readVariables('float-x','float-z');

figure, plot(x.',z.')