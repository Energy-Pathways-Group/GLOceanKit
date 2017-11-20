%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjustable parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wavelength = 1e3; % wavelength in meters
j = 1; % vertical mode number
epsilon = 0.2; % nonlinearity parameter
maxOscillations = 100; % Total number of oscillations, in periods
stratification = 'exponential';
z0 = [-10; -250; -625]; % initial particle positions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derived parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
latitude = 0.0;
g = 9.81;
k = 2*pi/wavelength;
[rho, N2, zIn] = InternalModes.StratificationProfileWithName(stratification);
zDomain = linspace(min(zIn),max(zIn),5000)';
im = InternalModes(rho,zIn,zDomain,latitude);
im.normalization = Normalization.uMax;
[F_out,G_out,h_out,omega_out] = im.ModesAtWavenumber(k);

F = F_out(:,j);
G = G_out(:,j);
h = h_out(j);
omega = omega_out(j);

U = epsilon*(omega/k);
period = abs(2*pi/omega);
t = (0:(1/40):maxOscillations)'*period;
fprintf('The wave period is set to %.1f hours.\n', period/3600)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Velocity field & numerical time stepping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = @(t, x) U*[cos(k*x(1)+omega*t)*interp1(zDomain,F,x(2)); k*h*sin(k*x(1) + omega*t)*interp1(zDomain,G,x(2))];

x = zeros(length(t),length(z0));
z = zeros(length(t),length(z0));
b = zeros(length(t),length(z0));

for i=1:length(z0)
    % Using ode113 for extremely high error tolerances
    options = odeset('RelTol',1e-5,'AbsTol',1e-5); % overkill
    [~, X] = ode113(u,t,[0 z0(i)],options);
    
    x(:,i) = X(:,1);
    z(:,i) = X(:,2);
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

