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

omega = wavemodel.f0;
alpha = 0;
% wavemodel.SetExternalWavesWithFrequencies(omega,alpha,j0,phi,U,'maxU');

wavemodel.InitializeWithPlaneWave(k0,l0,j0,U,sign);

% initial positions, center of the domain doesn't deal with boundaries
p0 = [Lx/2, Ly/2, -Lz/4];

% iniital position near zero, requires dealing with the boundary
p0 = [0, 0, -Lz/4];

timeStep = 15*60; % in seconds
maxTime = 2*pi/wavemodel.f0;
t_in = (0:timeStep:maxTime)';

% First let's do the adaptive time-stepping integrator
f = @(t,y) wavemodel.VelocityAtTimePositionVector(t,y);
[t,p] = ode23(f,t_in, p0);
x = p(:,1);
y = p(:,2);
z = p(:,3);

% Next, let's do fixed step size integrator.
cfl = 0.5;
deltaT = cfl*(wavemodel.x(2)-wavemodel.x(1))/U;
if (deltaT < t_in(2)-t_in(1))
   t_in = (0:60*ceil(deltaT/60):maxTime)';
end

% Try ode2, ode3, & ode4, depending on required accuracy.
p2 = ode2(f,t_in, p0');
x2 = p2(:,1);
y2 = p2(:,2);
z2 = p2(:,3);


figure
plot(x,y), hold on, plot(x2,y2)
legend('adaptive time step', 'fixed time step')