%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 256;
aspectRatio = 1/2;

Lx = 30e3;
Ly = aspectRatio*Lx;
Lz = 1300;

Nx = N;
Ny = aspectRatio*N;
Nz = N+1; % 2^n + 1 grid points, to match the Winters model, but 2^n ok too.

latitude = 25;
N0 = 5.2e-3/2; % Choose your stratification 7.6001e-04

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wvm = WaveVortexModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);

cg_x = wvm.cg_x;
cg_y = wvm.cg_y;
cg_z = wvm.cg_z;

a = squeeze(cg_z(:,1,:));
figure, pcolor(a), shading interp, colorbar('eastoutside')
b = squeeze(cg_x(:,1,:));
figure, pcolor(b), shading interp, colorbar('eastoutside')

return;

U = .2;
% boussinesq.InitializeWithPlaneWave(0,0,1,U,1); 
period = wvm.InitializeWithPlaneWave(10,0,1,U,1);  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the integrator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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