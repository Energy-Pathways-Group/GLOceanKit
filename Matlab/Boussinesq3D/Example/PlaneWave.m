%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 64;
aspectRatio = 1;

Lx = 100e3;
Ly = aspectRatio*100e3;
Lz = 1300;

Nx = N;
Ny = aspectRatio*N;
Nz = N+1; % 2^n + 1 grid points, to match the Winters model, but 2^n ok too.

latitude = 31;
N0 = 5.2e-3; % Choose your stratification 7.6001e-04

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

boussinesq = Boussinesq3DConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);

U = .2;
boussinesq.InitializeWithPlaneWave(2,2,1,U,1);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the integrator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = 2*pi/boussinesq.f0/100;
integrator = ArrayIntegrator(@(t,y0) boussinesq.NonlinearFluxAtTimeArray(t,y0),boussinesq.Y,dt);

nT=100;
totalEnergy = zeros(nT,1);
totalSpectralEnergy = zeros(nT,1);
% profile on
for i=1:100
   integrator.currentY = boussinesq.Y;
   boussinesq.Y = integrator.IncrementForward();
   totalEnergy(i) = boussinesq.totalEnergy;
   totalSpectralEnergy(i) = boussinesq.totalSpectralEnergy;
   if mod(i,10)==0
       boussinesq.summarizeEnergyContent();
   end
end
% profile viewer