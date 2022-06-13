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
model = WaveVortexModel(wvt);
model.nonlinearFlux = SingleModeQGPVE(model.wvt,shouldUseBeta=1);


% set initial positions for a bunch of floats
nTrajectories = 15000;
x = linspace(Lx/8,(7/8)*Lx, round((Lx/Ly)*sqrt(nTrajectories/(Lx/Ly))) );
y = linspace(Ly/8,(7/8)*Ly, round(sqrt(nTrajectories/(Lx/Ly))));
[xFloat,yFloat] = ndgrid(x,y);
xFloat = reshape(xFloat,1,[]);
yFloat = reshape(yFloat,1,[]);
nTrajectories = length(xFloat);
model.SetDrifterPositions(xFloat,yFloat,zeros(size(xFloat)),'qgpv');

[deltaT,advectiveDT,oscillatoryDT] = model.TimeStepForCFL(0.15);
finalTime = 75*86400;
nT = model.SetupIntegrator(advectiveDT, 86400, finalTime);

xFloatT = zeros(nT,nTrajectories);
yFloatT = zeros(nT,nTrajectories);
qgpvFloatT = zeros(nT,nTrajectories);
t = zeros(nT,1);
[xFloatT(1,:),yFloatT(1,:),~,qgpvFloatT(1,:)] = model.DrifterPositions;

while(model.t < finalTime)
    t(model.outputIndex) = model.integrateToNextOutputTime();
    [xFloatT(model.outputIndex,:),yFloatT(model.outputIndex,:),~,qgpvFloatT(model.outputIndex,:)] = model.DrifterPositions;
end

figure, pcolor(wvt.x,wvt.y,wvt.ssh.'), shading interp

return;

% model.linearDynamics = 1;


zFloat = linspace(-Lz,0,nTrajectories);

model.SetFloatPositions(xFloat,yFloat,zFloat,'rho_total');

% Set up the integrator
outputInterval = period/10;
deltaT = model.TimeStepForCFL(0.5,outputInterval);
finalTime = 3*period;
nT = model.SetupIntegrator(deltaT, outputInterval,finalTime);

% write the float trajectories to memory
xFloatT = zeros(nT,nTrajectories);
yFloatT = zeros(nT,nTrajectories);
zFloatT = zeros(nT,nTrajectories);
rhoFloatT = zeros(nT,nTrajectories);
t = zeros(nT,1);

[xFloatT(1,:),yFloatT(1,:),zFloatT(1,:),rhoFloatT(1,:)] = model.FloatPositions;

while(model.t < finalTime)
    t(model.outputIndex) = model.integrateToNextOutputTime();
    [xFloatT(model.outputIndex,:),yFloatT(model.outputIndex,:),zFloatT(model.outputIndex,:),rhoFloatT(model.outputIndex,:)] = model.FloatPositions;
end

figure, plot(xFloatT,zFloatT)
