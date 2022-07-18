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
model = WaveVortexModel(wvt,nonlinearFlux=SingleModeQGPVE(wvt,shouldUseBeta=1,u_damp=wvt.uMax));

% set initial positions for a bunch of floats
[xFloat,yFloat] = ndgrid(wvt.x(1:2:end),wvt.y(1:2:end));
xFloat = reshape(xFloat,1,[]);
yFloat = reshape(yFloat,1,[]);
nTrajectories = length(xFloat);
model.setDrifterPositions(xFloat,yFloat,[],'qgpv');

finalTime=150*86400;
nT = model.setupIntegrator(timeStepConstraint="advective", outputInterval=86400,finalTime=finalTime);

xFloatT = zeros(nT,nTrajectories);
yFloatT = zeros(nT,nTrajectories);
qgpvFloatT = zeros(nT,nTrajectories);
t = zeros(nT,1);
[xFloatT(1,:),yFloatT(1,:),~,tracked] = model.drifterPositions;
qgpvFloatT(model.outputIndex,:) = tracked.qgpv;

while(model.t < finalTime)
    t(model.outputIndex) = model.integrateToNextOutputTime();
    [xFloatT(model.outputIndex,:),yFloatT(model.outputIndex,:),~,tracked] = model.drifterPositions;
    qgpvFloatT(model.outputIndex,:) = tracked.qgpv;
end

figure
subplot(2,1,1)
pcolor(wvt.x/1000,wvt.y/1000,wvt.ssh.'), shading interp, axis equal
xlim([min(wvt.x) max(wvt.x)]/1000),ylim([min(wvt.y) max(wvt.y)]/1000)
colormap(lansey)
subplot(2,1,2)
plot(xFloatT/1000,yFloatT/1000), axis equal
xlim([min(wvt.x) max(wvt.x)]/1000),ylim([min(wvt.y) max(wvt.y)]/1000)
packfig(2,1,'rows')
% print('monopolefigure.png','-dpng','-r300')


