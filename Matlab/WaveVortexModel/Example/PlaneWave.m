%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 64;
aspectRatio = 1;

Lx = 50e3;
Ly = aspectRatio*Lx;
Lz = 1300;

Nx = N;
Ny = aspectRatio*N;
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



return

for iIncrement=1:totalIncrements
    self.ShowIntegrationTimeDiagnostics(iIncrement);

    self.IncrementForward();


end

% iTool.IntegrateToTime(wvm.inertialPeriod);

return

dt = period/50;
nT=5*50;
nTrajectories = 101;
totalEnergy = zeros(nT,1);
totalSpectralEnergy = zeros(nT,1);
totalEnergy(1) = wvm.totalEnergy;
totalSpectralEnergy(1) = wvm.totalSpectralEnergy;
x = zeros(nT,nTrajectories); y = zeros(nT,nTrajectories); z = zeros(nT,nTrajectories);
x(1,:) = Lx/2*ones(1,nTrajectories); y(1,:) = Ly/2*ones(1,nTrajectories); z(1,:) = linspace(-Lz,0,nTrajectories);

integrator = ArrayIntegrator(@(t,y0) wvm.NonlinearFluxWithParticlesAtTimeArray(t,y0),{wvm.Ap,wvm.Am,wvm.A0,x(1,:),y(1,:),z(1,:)},dt);

% profile on
for i=2:nT
%    integrator.currentY = wvm.Y;
   integrator.IncrementForward();
   wvm.Ap = integrator.currentY{1};
   wvm.Am = integrator.currentY{2};
   wvm.A0 = integrator.currentY{3};
   totalEnergy(i) = wvm.totalEnergy;
   totalSpectralEnergy(i) = wvm.totalSpectralEnergy;
   x(i,:) = integrator.currentY{4};
   y(i,:) = integrator.currentY{5};
   z(i,:) = integrator.currentY{6};
%    if mod(i,10)==0
       wvm.summarizeEnergyContent();
%    end
end
% profile viewer