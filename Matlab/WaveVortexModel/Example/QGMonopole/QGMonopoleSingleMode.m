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

wvt = WVTransformSingleMode([Lx, Ly], [Nx, Ny], h=0.8, latitude=latitude);

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
model = WVModel(wvt,nonlinearFlux=SingleModeQGPVE(wvt,shouldUseBeta=1,u_damp=wvt.uMax));

if model.nonlinearFluxOperation.beta > 0
    beta = 2 * 7.2921E-5 * cos( wvt.latitude*pi/180. ) / 6.371e6;
    outputVar = WVVariableAnnotation('qgpv',{'x','y','z'},'1/s', 'quasigeostrophic potential vorticity');
    f = @(wvt) -wvt.transformToSpatialDomainWithF( (wvt.Omega .* wvt.Omega ./ (wvt.h * wvt.f)) .*wvt.A0t) + beta*wvt.Y;
    wvt.addOperation(WVOperation('qgpv',outputVar,f),overwriteExisting=1);
    f = @(wvt) wvt.diffX(wvt.v) - wvt.diffY(wvt.u) - (wvt.f)^2/(wvt.h * wvt.g)*wvt.psi + beta*wvt.Y;
end


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

%%
[~,index] = max(std(qgpvFloatT,0,1));

figure
tiledlayout(2,1)
nexttile
plot(t/86400,qgpvFloatT(:,index)/wvt.f);
ylabel('qgpv (f_0)')
xlabel('time (days)')
nexttile
plot(xFloatT(:,index)/1000,yFloatT(:,index)/1000)


