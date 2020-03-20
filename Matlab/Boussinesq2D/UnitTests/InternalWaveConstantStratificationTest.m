%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Problem parameters

epsilon = 0.15;
nPeriods = 10;
zeta0 = 100;

N0 = 5.2e-3;
D = 1300;
k = 2*pi/100;
j = 1;
m = pi*j/D;
omega = N0*k/sqrt(k*k+m*m);
rho0 = 1025;
g = 9.81;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c_p = omega/k;
U = c_p*epsilon;

t = linspace(0,nPeriods*(2*pi/omega),5000).';

fprintf('Fluid velocity set to %.2f cm/s\n',U*100);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Eulerian velocity

norm = -1;
psi = @(t,x,z) norm*(U/m)*cos(k*x+omega*t).*sin(m*z);
u = @(t,x,z) norm*U*cos(k*x+omega*t).*cos(m*z);
w = @(t,x,z) norm*U*(k/m)*sin(k*x+omega*t).*sin(m*z);
b = @(t,x,z) norm*N0*N0*(U*k/(m*omega))*cos(k*x+omega*t).*sin(m*z);
uw = @(t,x) [u(t,x(1),x(2));w(t,x(1),x(2))];
% uw = @(t,x) [U*cos(k*x(1)+omega*t)*cos(m*x(2));U*(k/m)*sin(k*x(1)+omega*t)*sin(m*x(2))];

% leaving off rho0 - ...
% rho = @(t,x,z) (rho0/g)*N0*N0*(z+(U*k/(m*omega))*cos(k*x+omega*t).*sin(m*z));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Numerical Lagrangian velocity

options = odeset('RelTol',1e-5,'AbsTol',1e-5); % overkill
sol = ode113(uw,[min(t) max(t)],[0,zeta0],options);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Setup the 2D Bousinesq model with the same parameters

Lx = 4*2*pi/k;
Lz = D;
Nx = 256;
Nz = 65;
N2 = N0*N0;

model = Boussinesq2DConstantStratification([Lx, Lz], [Nx, Nz], linspace(-Lz,0,Nz).', N2);
[omega,h] = model.InitializeWithPlaneWave(k,j,U);
model.setParticlePositions(Lx/2,zeta0-Lz);

u_model=model.u;
w_model=model.w;
b_model=model.b;

i=4;
x_i = model.x(i);
s = model.z + Lz;


% figure
% subplot(1,3,1)
% plot(u_model(i,:),model.z), hold on
% plot(u(0,x_i,s),model.z)
% subplot(1,3,2)
% plot(w_model(i,:),model.z), hold on
% plot(w(0,x_i,s),model.z)
% subplot(1,3,3)
% plot(b_model(i,:),model.z), hold on
% plot(b(0,x_i,s),model.z)

rel_error = @(model,truth) max(abs(model(:)-truth(:)))./max(abs(truth(:)));

n = 100;
t = zeros(n,1);
xi = zeros(n,1);
zeta = zeros(n,1);
u_err = zeros(n,1);
w_err = zeros(n,1);
b_err = zeros(n,1);
for i=1:n
    model.StepForwardToTime(i*model.dt);
    t(i) = model.t;
    xi(i,:) = model.xi;
    zeta(i,:) = model.zeta;
    
    u_err(i)=rel_error(model.u,u(model.t,model.X,model.Z+Lz));
    w_err(i)=rel_error(model.w,w(model.t,model.X,model.Z+Lz));
    b_err(i)=rel_error(model.b,b(model.t,model.X,model.Z+Lz));
end

figure
plot(t,[u_err,w_err,b_err])

p = deval(t,sol);

x = p(1,:).';
z = p(2,:).';

figure
plot(x-x(1),z-Lz), hold on
plot(xi-xi(1),zeta)