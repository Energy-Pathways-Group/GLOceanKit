Lx = 15e3;
Ly = 15e3;
Lz = 5000;

Nx = 64;
Ny = 64;
Nz = 64;

latitude = 31;
N0 = 5.2e-3/2; % Choose your stratification 7.6001e-04

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wavemodel = InternalWaveModelSlow([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);

k0 = 2; % k=0..Nx/2
l0 = 0; % l=0..Ny/2
j0 = 2; % j=1..nModes, where 1 indicates the 1st baroclinic mode
U = 0.01; % m/s
sign = 1;

wavemodel.InitializeWithPlaneWave(k0,l0,j0,U,sign);

timeStep = 1; % in seconds
maxTime = 100;
t = (0:timeStep:maxTime)';

tic
for iTime=1:length(t)
    [u,v]=wavemodel.VelocityFieldAtTime(t(iTime));
    [w,zeta] = wavemodel.VerticalFieldsAtTime(t(iTime));
end
toc