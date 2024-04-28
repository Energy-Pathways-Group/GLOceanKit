N = 128;

Lx = 750e3;
Ly = 750e3;
Lz = 4000;

Nx = N;
Ny = N;
nModes = 40;

latitude = 31;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N0 = 5.2e-3; % reference buoyancy frequency, radians/seconds
rho0 = 1025; g = 9.81;
L_gm = 1.3e3; % thermocline exponential scale, meters
rhoFunc = @(z) rho0*(1 + L_gm*N0*N0/(2*g)*(1 - exp(2*z/L_gm)));
% N2Func = @(z) N0*N0*exp(2*z/L_gm);
zIn = [-2000 0];

wvm = WaveVortexModelHydrostatic([Lx, Ly, max(zIn)-min(zIn)], [Nx, Ny, nModes], latitude, rhoFunc);

x0 = (max(wvm.x)-min(wvm.x))/2;
y0 = (max(wvm.y)-min(wvm.y))/2;

A = -1.32e-5; % s^{-1}
alpha = 8e-10; % m^{-2}
beta = 8.2e-6; % m^{-2}
gamma = 0.01; % m^{-1}
psi = @(x,y,z) -(A/(2*alpha))*exp(-alpha*((x-x0).*(x-x0)+(y-y0).*(y-y0))-beta*z.*z);

wvm.setGeostrophicStreamfunction(psi);
[u,v,eta] = wvm.VariableFieldsAtTime(0, 'u','v','eta');

x = wvm.x;
y = wvm.y;
z = wvm.z;
f = wvm.f;

[X,Y,Z] = ndgrid(wvm.x,wvm.y,wvm.z);
psiG = psi(X,Y,Z);
psiGy = -DiffFourier(y,psiG,1,2);
figure, pcolor(x/1e3,z/1e3,squeeze(psiGy(:,Ny/2,:)).'), shading interp; colorbar('eastoutside')
figure, pcolor(x/1e3,y/1e3,squeeze(psiGy(:,:,35)).'), shading interp; colorbar('eastoutside')
max(abs(psiGy(:)))


max(abs(u(:)))
figure, pcolor(x/1e3,z/1e3,squeeze(u(:,Ny/2,:)).'), shading interp; colorbar('eastoutside')
figure, pcolor(x/1e3,y/1e3,squeeze(u(:,:,35)).'), shading interp; colorbar('eastoutside')

zeta_z = (DiffFourier(x,v,1,1) - DiffFourier(y,u,1,2))/f;

figure, pcolor(x/1e3,z/1e3,squeeze(zeta_z(:,Ny/2,:)).'), shading interp; colorbar('eastoutside')
figure, pcolor(x/1e3,y/1e3,squeeze(zeta_z(:,:,35)).'), shading interp; colorbar('eastoutside')


% B = 0.2; % m/s
