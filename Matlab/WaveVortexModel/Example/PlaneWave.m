%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 64;
aspectRatio = 1;

Lx = 500e3;
Ly = aspectRatio*Lx;
Lz = 1300;

Nx = N;
Ny = aspectRatio*N;
Nz = N+1; % 2^n + 1 grid points, to match the Winters model, but 2^n ok too.

latitude = 25;
N0 = 5.2e-3; % Choose your stratification 7.6001e-04

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wvm = WaveVortexModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);

U = .2;
% boussinesq.InitializeWithPlaneWave(0,0,1,U,1); 
wvm.InitializeWithPlaneWave(2,2,1,U,1);  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the integrator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = 2*pi/wvm.f0/50;
integrator = ArrayIntegrator(@(t,y0) wvm.NonlinearFluxAtTimeArray(t,y0),{wvm.Ap,wvm.Am,wvm.A0},dt);

nT=20*50;
totalEnergy = zeros(nT,1);
totalSpectralEnergy = zeros(nT,1);
% profile on
for i=1:nT
%    integrator.currentY = wvm.Y;
   integrator.IncrementForward();
   wvm.Ap = integrator.currentY{1};
   wvm.Am = integrator.currentY{2};
   wvm.A0 = integrator.currentY{3};
   totalEnergy(i) = wvm.totalEnergy;
   totalSpectralEnergy(i) = wvm.totalSpectralEnergy;
%    if mod(i,10)==0
       wvm.summarizeEnergyContent();
%    end
end
% profile viewer