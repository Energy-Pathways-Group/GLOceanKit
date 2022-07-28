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

wvt = WVTransformConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], N0,latitude=latitude);

U = .2; phi=0;
omega = wvt.initWithWaveModes(10,0,1,phi,U,1);
period = 2*pi/omega;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the integrator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize the integrator with the model
%model = WVModel(wvt,nonlinearFlux=NonlinearBoussinesqWithReducedInteractionMasks(wvt));
model = WVModel(wvt);

% set initial positions for a bunch of floats
nTrajectories = 101;
xFloat = Lx/2*ones(1,nTrajectories);
yFloat = Ly/2*ones(1,nTrajectories);
zFloat = linspace(-Lz,0,nTrajectories);

model.setFloatPositions(xFloat,yFloat,zFloat,'rho_total');

% Set up the integrator
nT = model.setupIntegrator(timeStepConstraint="oscillatory", outputInterval=period/10,finalTime=3*period);

% write the float trajectories to memory
xFloatT = zeros(nT,nTrajectories);
yFloatT = zeros(nT,nTrajectories);
zFloatT = zeros(nT,nTrajectories);
rhoFloatT = zeros(nT,nTrajectories);
t = zeros(nT,1);

[xFloatT(1,:),yFloatT(1,:),zFloatT(1,:),tracked] = model.floatPositions;
rhoFloatT(1,:) = tracked.rho_total;

while(model.outputIndex < nT)
    t(model.outputIndex) = model.integrateToNextOutputTime();
    [xFloatT(model.outputIndex,:),yFloatT(model.outputIndex,:),zFloatT(model.outputIndex,:),tracked] = model.floatPositions;
    rhoFloatT(model.outputIndex,:) = tracked.rho_total;
end

figure, plot(xFloatT,zFloatT)
