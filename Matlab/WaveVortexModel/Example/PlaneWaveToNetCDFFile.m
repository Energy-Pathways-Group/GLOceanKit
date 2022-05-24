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

wvt = WaveVortexTransformConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);

U = .2;
period = wvt.InitializeWithPlaneWave(10,0,1,U,1);  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the integrator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize the integrator with the model
model = WaveVortexModel(wvt);

% set initial positions for a bunch of floats
nTrajectories = 101;
xFloat = Lx/2*ones(1,nTrajectories);
yFloat = Ly/2*ones(1,nTrajectories);
zFloat = linspace(-Lz,0,nTrajectories);
model.SetFloatPositions(xFloat,yFloat,zFloat);

% Set up the integrator
outputInterval = period/10;
deltaT = model.TimeStepForCFL(0.5,outputInterval);
finalTime = 3*period;
nT = model.SetupIntegrator(deltaT, outputInterval,finalTime);

model.CreateNetCDFFileForModelOutput('PlaneWaveWithFloats.nc',shouldOverwriteExisting=1);

model.IntegrateToTime(finalTime);

ncfile = model.ncfile;
% [x,y,z] = ncfile.FloatPositions();
[x,z] = ncfile.readVariables('float-x','float-z');

figure, plot(x.',z.')