Lx = 15e3;
Ly = 15e3;
Lz = 5000;

Nx = 64;
Ny = 64;
Nz = 64;

latitude = 31;
N0 = 5.2e-3; % Choose your stratification

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wavemodel = InternalWaveModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);

j0 = 1; % j=1..nModes, where 1 indicates the 1st baroclinic mode
U = 0.1; % m/s
sign = 1;
phi = 0;
k0 = 0;
l0 = 0;

omega = 4*wavemodel.f0;
alpha = 0;
wavemodel.SetExternalWavesWithFrequencies(omega,alpha,j0,phi,U,'maxU');

% wavemodel.InitializeWithPlaneWave(k0,l0,j0,U,sign);

% initial positions
p0 = [Lx/2, Ly/2, -Lz/4];

timeStep = 15*60; % in seconds
maxTime = 2*pi/wavemodel.f0;

[t,p] = ode23(@(t,y) wavemodel.VelocityAtTimePositionVector(t,y),[0 maxTime], p0);
x = p(:,1);
y = p(:,2);
z = p(:,3);

figure
plot(x,y)