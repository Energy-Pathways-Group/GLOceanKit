Lx = 10e3;
Lz = 1300;
Nx = 256;
Nz = 129;
z = linspace(-Lz,0,Nz).';
N2 = (5.2e-3)^2;

model = Boussinesq2D([Lx, Lz], [Nx, Nz], z, N2);

U = 0.2;
j = 1;
k = model.k(4);
[omega,h] = model.InitializeWithPlaneWave(k,j,U);
epsilon = U/(omega/k)

figure
subplot(1,2,1)
pcolor(model.X,model.Z,model.psi_n), shading interp
subplot(1,2,2)
pcolor(model.X,model.Z,model.b_n), shading interp

figure
pcolor(model.X,model.Z,model.nabla2_psi_n), shading interp

cfl = 0.25;


f = @(t,y0) model.linearFluxCat(y0);
integrator = Integrator(f,cat(3,model.nabla2_psi_n,model.b),dt);