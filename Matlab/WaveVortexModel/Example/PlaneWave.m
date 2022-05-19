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

wvm = WaveVortexModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);

U = .2;
period = wvm.InitializeWithPlaneWave(10,0,1,U,1);  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the integrator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize the integrator with the model
iTool = WaveVortexModelIntegrationTools(wvm);

% set initial positions for a bunch of floats
nTrajectories = 101;
iTool.xFloat = Lx/2*ones(1,nTrajectories);
iTool.yFloat = Ly/2*ones(1,nTrajectories);
iTool.zFloat = linspace(-Lz,0,nTrajectories);

% Set up the integrator
outputInterval = period/10;
deltaT = iTool.TimeStepForCFL(0.5,outputInterval);
finalTime = 3*period;
nT = iTool.SetupIntegrator(deltaT, outputInterval,finalTime);

% write the float trajectories to memory
xFloatT = zeros(nT,nTrajectories);
yFloatT = zeros(nT,nTrajectories);
zFloatT = zeros(nT,nTrajectories);
t = zeros(nT,1);

xFloatT(1,:) = iTool.xFloat;
yFloatT(1,:) = iTool.yFloat;
zFloatT(1,:) = iTool.zFloat;

while(iTool.t < finalTime)
    t(iTool.outputIndex) = iTool.integrateToNextOutputTime();
    xFloatT(iTool.outputIndex,:) = iTool.xFloat;
    yFloatT(iTool.outputIndex,:) = iTool.yFloat;
    zFloatT(iTool.outputIndex,:) = iTool.zFloat;
end

figure, plot(xFloatT,zFloatT)
