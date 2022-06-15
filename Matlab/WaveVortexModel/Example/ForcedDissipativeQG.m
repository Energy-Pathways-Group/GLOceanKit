%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lx = 250e3;
Ly = 250e3;

Nx = 128;
Ny = 128;

latitude = 25;

wvt = WaveVortexTransformSingleMode([Lx, Ly], [Nx, Ny], h=0.8, latitude=latitude);

dk = (wvt.k(2)-wvt.k(1));
k_f = 20*dk;
f_zeta = 0.01*(wvt.f0)^2;
k_r = 4*dk;

epsilon = 10*(f_zeta.^(3/2))/(k_f.*k_f); % m^2/s^3
r = 0.04*(epsilon*k_r^2)^(1/3); % 1/s
u_rms = (epsilon/k_r)^(1/3); % m/s
nu = (3/2)*(wvt.x(2)-wvt.x(1))*u_rms; % m^2/s

fdFlux = SingleModeForcedDissipativeQGPVE(wvt,k_f=k_f,f_zeta=f_zeta,r=r,nu=nu);
model = WaveVortexModel(wvt,nonlinearFlux=fdFlux);

deltaT = (wvt.x(2)-wvt.x(1))*0.25/u_rms;
model.setupIntegrator(deltaT=deltaT, outputInterval=86400);

% model.integrateOneTimeStep();

% model.integrateToNextOutputTime();
model.integrateToTime(1000*86400);

figure, pcolor(wvt.x,wvt.y,wvt.ssh.'), shading interp

A02 = wvt.A0_TE_factor .* (wvt.A0.*conj(wvt.A0));
GeostrophicEnergyK = wvt.transformToRadialWavenumber(A02);
figure
plot(wvt.kRadial,GeostrophicEnergyK), xlog, ylog
vlines(k_f,'g--')

figure, pcolor(log10(abs(fftshift(fdFlux.damp)))), shading interp
figure, pcolor(log10(abs(fftshift(fdFlux.F)))), shading interp

DampK = wvt.transformToRadialWavenumber(fdFlux.damp);