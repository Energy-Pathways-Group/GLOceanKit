%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjustable parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wavelength = 1000; % wavelength in meters
j = 1; % vertical model number
epsilon = 0.1; % nonlinearity parameter
maxOscillations = 100; % Total number of oscillations, in periods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fixed parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = 500;
N2 = (5.3e-3)^2;
rho0 = 1025;
g = 9.81;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derived parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 2*pi/wavelength;
m = j*pi/D;
h = N2/(k*k + m*m)/g;
omega = sqrt(g*h*k*k);
U = epsilon*(omega/k);
drhodz = -(rho0/g)*N2;

period = abs(2*pi/omega);
t = (0:(1/40):maxOscillations)'*period;
fprintf('The wave period is set to %.1f hours.\n', period/3600)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Velocity field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = @(t,x) U*[cos(k*x(1)+omega*t)*cos(m*(x(2)+D));  (k/m)*sin(k*x(1)+omega*t)*sin(m*(x(2)+D))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% density
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho = @(t,x) -(rho0*N2/g)*(x(:,2) + (epsilon/m)*cos(k*x(:,1)+omega*t).*sin(m*(x(:,2)+D)) );

depths = linspace(-D+5,-5,5)';
x = zeros(length(t),length(depths));
z = zeros(length(t),length(depths));
b = zeros(length(t),length(depths));

for i=1:length(depths)
    % Using ode113 for extremely high error tolerances
    options = odeset('RelTol',1e-12,'AbsTol',1e-12); % overkill
    [T, X] = ode113(u,t,[0 depths(i)],options);
    
    x(:,i) = X(:,1);
    z(:,i) = X(:,2);
    b(:,i) = rho(t,[x(:,i),z(:,i)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Theoretical 2nd order Stokes drift
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U_stokes = @(z) -(U*U*k/(2*omega))*(cos(m*(z-D)).*cos(m*(z-D)) - 1/(h*m*m*9.81)*(N2-omega*omega)*sin(m*(z-D)).*sin(m*(z-D)) );

figure
plot(x,z), hold on
z_stokes = linspace(-D,0,500);
plot(U_stokes(z_stokes)*max(t),z_stokes)

return

u_vec = @(t,x) U*[cos(k*x(:,1)+omega*t).*cos(m*(x(:,2)+D)),  (k/m)*sin(k*x(:,1)+omega*t).*sin(m*(x(:,2)+D))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare Lagragian velocity to Eulerian velocity at the particle position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i = 3;
dt = t(2)-t(1);
xi = x(:,i);
zeta = z(:,i);
den = rho0 + b(:,i);

rho_i = rho0 - rho(0,[0,depths(i)]);
zeta_bar = (g/N2)*(rho_i - rho0)/rho0;
figure
plot(t,epsilon*epsilon*(sin(m*zeta).^2) - m*m*(zeta_bar - zeta).^2)

xi_t = vdiff(dt,xi,1);
zeta_t = vdiff(dt,zeta,1);

figure, plot(t, xi_t.^2 + zeta_t.^2 )

u_particle = u_vec(t,[xi,zeta]);
figure
subplot(2,1,1)
plot(t,u_particle(:,1)), hold on
plot(t,xi_t)
subplot(2,1,2)
plot(t,u_particle(:,2)), hold on
plot(t,zeta_t)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare Lagragian velocity to Eulerian velocity at the particle position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zeta_tt = vdiff(dt,zeta_t,1);
figure, plot(t,zeta_tt), hold on
plot(t,(U*U*k*k/m)*sin(m*(zeta+D)).*cos(m*(zeta+D)) + U*omega*cos(k*xi+omega*t).*sin(m*(zeta+D)))

return

figure, plot(t,sin(m*z)./sin(m*z(1,:)))

dt = t(2)-t(1);
zeta = z(:,5);
zeta_t = vdiff(dt,zeta,1);



% plot(t,(U*U*k*k/(2*m))*sin(2*m*zeta))

a = ((U*k/m)^2)*(sin(m*zeta).^2 - sin(m*(zeta(1)))^2);
b = (diff(zeta)./diff(t)).^(2);
figure, plot(a), hold on, plot(b)

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analytical solution?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i = 2;
[sn,cn,dn] = ellipj(-(U*k)*t,sin(m*depths(i))^2);
zeta = asin(sin(m*depths(i))*1./dn)/m;
figure, plot(t,zeta), hold on, plot(t,z(:,i))
figure, plot(t, sin(m*z(:,i))), hold on, plot(t, sin(m*depths(i))*cn./dn)