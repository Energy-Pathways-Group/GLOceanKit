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
% Setup the model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N0 = 5.2e-3; % reference buoyancy frequency, radians/seconds
rho0 = 1025; g = 9.81;
L_gm = 1.3e3; % thermocline exponential scale, meters
rhoFunc = @(z) rho0*(1 + L_gm*N0*N0/(2*g)*(1 - exp(2*z/L_gm)));
zIn = [-2000 0];

wvm = WaveVortexModelHydrostatic([Lx, Ly, max(zIn)-min(zIn)], [Nx, Ny, nModes], latitude, rhoFunc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Add a geostrophic stream function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x0 = (max(wvm.x)-min(wvm.x))/2;
y0 = (max(wvm.y)-min(wvm.y))/2;

A = -1.32e-5; % s^{-1}
alpha = 8e-10; % m^{-2}
beta = 8.2e-6; % m^{-2}
psi = @(x,y,z) -(A/(2*alpha))*exp(-alpha*((x-x0).*(x-x0)+(y-y0).*(y-y0))-beta*z.*z);

wvm.SetGeostrophicStreamfunction(psi);
[u,v,eta] = wvm.VariableFieldsAtTime(0, 'u','v','eta');

x = wvm.x;
y = wvm.y;
z = wvm.z;
f0 = wvm.f0;
% zeta_z = (DiffFourier(x,v,1,1) - DiffFourier(y,u,1,2))/f0;
% figure, pcolor(x/1e3,z/1e3,squeeze(zeta_z(:,Ny/2,:)).'), shading interp; colorbar('eastoutside')
% figure, pcolor(x/1e3,y/1e3,squeeze(zeta_z(:,:,35)).'), shading interp; colorbar('eastoutside')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Add a an inertial function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B = 0.2; % m/s
gamma = 0.01; % m^{-1}
u_IO = @(z) B*exp(gamma*z);
v_IO = @(z) zeros(size(z));

wvm.SetInertialMotions(u_IO,v_IO);

% N=5;tmp=cat(1,squeeze(Ubar(1,1,1:end-N)),zeros(N,1));figure, plot(self.PFinv * tmp,self.z)

% [u,v,eta] = wvm.VariableFieldsAtTime(0, 'u','v','eta');
% figure, pcolor(x/1e3,z/1e3,squeeze(u(:,Ny/2,:)).'), shading interp; colorbar('eastoutside')
% figure, pcolor(x/1e3,y/1e3,squeeze(u(:,:,floor(0.9*nModes))).'), shading interp; colorbar('eastoutside')

wvm.summarizeEnergyContent;

dt = 2*pi/wvm.f0/50;
integrator = ArrayIntegrator(@(t,y0) wvm.NonlinearFluxAtTimeArray(t,y0),{wvm.Ap,wvm.Am,wvm.A0},dt);

nT=20*50;
totalEnergy = zeros(nT,1);
totalSpectralEnergy = zeros(nT,1);
% profile on
for i=2:nT
   integrator.IncrementForward();
   wvm.Ap = integrator.currentY{1};
   wvm.Am = integrator.currentY{2};
   wvm.A0 = integrator.currentY{3};
   totalEnergy(i) = wvm.totalEnergy;
   totalSpectralEnergy(i) = wvm.totalSpectralEnergy;
    if mod(i,10)==0
       wvm.summarizeEnergyContent();
    end
end

[k,j,ApKJ,AmKJ] = wvm.ConvertToWavenumberAndMode(abs(wvm.Ap).^2,abs(wvm.Am).^2);
figure, plot(k,sum(ApKJ,2)), hold on, plot(k,sum(AmKJ,2)), ylog