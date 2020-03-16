Lx = 10e3;
Lz = 1300;
Nx = 256;
Nz = 65;
z = linspace(-Lz,0,Nz).';
N2 = (5.2e-3)^2;

model = Boussinesq2D([Lx, Lz], [Nx, Nz], z, N2);

U = 0.2;
j = 1;
k = model.k(4);
m = j*pi/Lz;
[omega,h] = model.InitializeWithPlaneWave(k,j,U);
epsilon = U/(omega/k)
W = (k/m)*U;

nParticles = 5;
model.setParticlePositions(Lx*ones(1,nParticles)/2,(Lz/6)*(1:nParticles) - Lz)

% figure
% subplot(1,2,1)
% pcolor(model.X,model.Z,model.psi_n), shading interp
% subplot(1,2,2)
% pcolor(model.X,model.Z,model.b_n), shading interp
% 
% figure
% pcolor(model.X,model.Z,model.nabla2_psi_n), shading interp

% cfl = 0.25;
% 
% dt_u = 0.25*(model.x(2)-model.x(1))/U;
% dt_w = 0.25*(model.z(2)-model.z(1))/W;
% 
% fprintf('cfl condition for (u,w)=>dt=(%.1f,%.1f) seconds for period of %.1f seconds\n',dt_u,dt_w,2*pi/omega);
% 
% dt = min(dt_u,dt_w);
% n = round(2*pi/omega/dt);
% dt = 2*pi/omega/n;

% f = @(t,y0) model.linearFluxCat(y0);
% integrator = Integrator(f,cat(2,model.nabla2_psi_n,model.b_n),dt);
% y = integrator.StepForwardToTime(30*dt);
% nabla2_psi = y(:,1:Nz);
% b = y(:,(Nz+1):end);

% f = @(t,y0) model.linearFlux(y0);
% y0 = {model.nabla2_psi_n;model.b_n};
% integrator = ArrayIntegrator(f,y0,dt);
% y = integrator.StepForwardToTime(30*dt);
% nabla2_psi = y{1};
% b = y{2};

% b0 = model.b;
% model.StepForwardToTime(30*dt);
% bn = model.b;
% 
% figure
% subplot(1,2,1)
% pcolor(model.X,model.Z,b0), shading interp
% subplot(1,2,2)
% pcolor(model.X,model.Z,bn), shading interp

figure
pcolor(model.X,model.Z,model.u), shading interp

n = 100;
maxU = zeros(n,1);
maxW = zeros(n,1);
xi = zeros(n,nParticles);
zeta = zeros(n,nParticles);
for i=1:n
    model.StepForwardToTime(i*model.dt);
    [u,w] = model.uw;
    maxU(i) = max(u(:));
    maxW(i) = max(w(:));
    xi(i,:) = model.xi;
    zeta(i,:) = model.zeta;
end

figure
plot([maxU, maxW])

figure
plot(xi,zeta)