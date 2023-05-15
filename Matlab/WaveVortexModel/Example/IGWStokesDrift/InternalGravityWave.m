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
% model = WVModel(wvt,nonlinearFlux=NonlinearBoussinesqWithReducedInteractionMasks(wvt));
model = WVModel(wvt);
% model = WVModel(wvt,nonlinearFlux=WaveWaveConstantN(wvt));

% set initial positions for a bunch of floats
nTrajectories = 101;
xFloat = wvt.Lx/2*ones(1,nTrajectories);
yFloat = wvt.Ly/2*ones(1,nTrajectories);
zFloat = linspace(-wvt.Lz,0,nTrajectories);

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
