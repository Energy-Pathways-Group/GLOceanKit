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


outputVar = WVVariableAnnotation('zeta',{'x','y','z'},'1/s^2', 'vertical component of relative vorticity');
wvt.addOperation(WVOperation('zeta',outputVar,@(wvt) wvt.diffX(wvt.v) - wvt.diffY(wvt.u)));

outputVar = WVVariableAnnotation('sigma_n',{'x','y','z'},'1/s^2', 'normal strain');
wvt.addOperation(WVOperation('sigma_n',outputVar,@(wvt) wvt.diffX(wvt.u) - wvt.diffY(wvt.v)));

outputVar = WVVariableAnnotation('sigma_s',{'x','y','z'},'1/s^2', 'shear strain');
wvt.addOperation(WVOperation('sigma_s',outputVar,@(wvt) wvt.diffX(wvt.v) + wvt.diffY(wvt.u)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the integrator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize the integrator with the model
model = WVModel(wvt,nonlinearFlux=SingleModeQGPVE(wvt,shouldUseBeta=1,u_damp=wvt.uMax));

% set initial positions for a bunch of floats
[xFloat,yFloat] = ndgrid(wvt.x(1:2:end),wvt.y(1:2:end));
xFloat = reshape(xFloat,1,[]);
yFloat = reshape(yFloat,1,[]);
nTrajectories = length(xFloat);
model.setDrifterPositions(xFloat,yFloat,[],'ssh','qgpv','u','v','zeta','sigma_n','sigma_s');

model.setupIntegrator(timeStepConstraint="advective", outputInterval=86400);

model.createNetCDFFileForModelOutput('QGMonopole.nc',shouldOverwriteExisting=1);
model.setNetCDFOutputVariables('A0','ssh','qgpv','u','v','zeta','sigma_n','sigma_s');
model.integrateToTime(1*86400);

% ncfile = model.ncfile;
% [x,y] = ncfile.readVariables('drifter-x','drifter-y');
% qgpv = ncfile.readVariables('drifter-qgpv');
% 
% figure, plot(x.',y.')
