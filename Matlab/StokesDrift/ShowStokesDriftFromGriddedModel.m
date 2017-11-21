%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjustable parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wavelength = 1e3; % wavelength in meters
j = 1; % vertical mode number
epsilon = 0.2; % nonlinearity parameter
maxOscillations = 5; % Total number of oscillations, in periods
stratification = 'exponential';
z0 = [-10; -250; -625]; % initial particle positions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derived parameters & model setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rho, N2, zIn] = InternalModes.StratificationProfileWithName(stratification);
latitude = 0.0;
k = 2*pi/wavelength;
wavenumber = 4;
Nx = 128;
Ny = 2;
Nz = 513; % Must include end point to advect at the surface, so use 2^N + 1
Lx = wavenumber*wavelength;
Ly = Ny*Lx/Nx;
Lz = max(zIn)-min(zIn);
zDomain = linspace(min(zIn),max(zIn),Nz)';

wavemodel = InternalWaveModelArbitraryStratification([Lx, Ly, Lz], [Nx, Ny, Nz], rho, zDomain, j+1, latitude);

% Pull out the wave modes to sketch out the theory for later (and get the
% wave period).
wavemodel.internalModes.normalization = Normalization.uMax;
[F_out,G_out,h_out,omega_out] = wavemodel.internalModes.ModesAtWavenumber(k);
F = F_out(:,j);
G = G_out(:,j);
h = h_out(j);
omega = omega_out(j);
U = epsilon*(omega/k);

period = wavemodel.InitializeWithPlaneWave(wavenumber,0,j,U,1);

t = (0:(1/40):maxOscillations)'*period;
fprintf('The wave period is set to %.1f hours.\n', period/3600)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Velocity field & numerical time stepping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = @(t, x) wavemodel.VelocityAtTimePositionVector(t,x);
x = zeros(length(t),length(z0));
y = zeros(length(t),length(z0));
z = zeros(length(t),length(z0));
b = zeros(length(t),length(z0));

for i=1:length(z0)
    % Using ode113 for extremely high error tolerances
    options = odeset('RelTol',1e-5,'AbsTol',1e-5); % overkill
    [~, X] = ode113(u,t,[0 0 z0(i)],options);
    
    x(:,i) = X(:,1);
    y(:,i) = X(:,2);
    z(:,i) = X(:,3);
%     b(:,i) = rho(t,[x(:,i),z(:,i)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Theoretical Stokes drift
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zTheory = linspace(-1000,0,1000)';
U1 = (U*U*k)/(2*omega)*(interp1(zDomain,F,zTheory)).^2;
U2 = -(U*U*k)/(2*omega)*(h/g*(N2(zTheory)-omega*omega)).*(interp1(zDomain,G,zTheory)).^2;

figure
plot(-(U1+U2)*max(t),zTheory)
hold on
plot(x,z,'LineWidth',2)

