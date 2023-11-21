%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lx = 3*2000e3;
Ly = 3*1000e3;

Nx = 2*256;
Ny = 2*128;

latitude = 25;

wvt = WVTransformSingleMode([Lx, Ly], [Nx, Ny], h=0.8, latitude=latitude);

x0 = 3*Lx/4;
y0 = Ly/2;
A = 0.15;
L = 80e3;
wvt.setSSH(@(x,y) A*exp( - ((x-x0).^2 + (y-y0).^2)/L^2) );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the integrator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize the integrator with the model
model = WVModel(wvt,nonlinearFlux=QGPVE(wvt,shouldUseBeta=1,u_damp=wvt.uMax));

model.setupIntegrator(timeStepConstraint="advective");
model.integrateToTime(3*365*86400);

%%
model.integrateToTime(4*365*86400);

figure, pcolor(wvt.x/1000,wvt.y/1000,wvt.ssh.'), shading interp, axis equal
colorbar('eastoutside')
xlim([min(wvt.x) max(wvt.x)]/1000); ylim([min(wvt.y) max(wvt.y)]/1000);
