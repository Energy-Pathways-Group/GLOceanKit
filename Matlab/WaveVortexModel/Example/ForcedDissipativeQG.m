%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lx = 50e3;
Ly = 50e3;

Nx = 128;
Ny = 128;

latitude = 25;

wvt = WaveVortexTransformSingleMode([Lx, Ly], [Nx, Ny], h=0.8, latitude=latitude);

dk = (wvt.k(2)-wvt.k(1));
k_f = 20*dk;
k_r = 1*dk;

% f_zeta = 0.01*(wvt.f0)^2;
% epsilon = 10*(f_zeta.^(3/2))/(k_f.*k_f); % m^2/s^3
% r = 0.04*(epsilon*k_r^2)^(1/3); % 1/s
% u_rms = (epsilon/k_r)^(1/3); % m/s
% nu = (3/2)*(wvt.x(2)-wvt.x(1))*u_rms; % m^2/s
% 
% fdFlux = SingleModeForcedDissipativeQGPVE(wvt,k_f=k_f,f_zeta=f_zeta,r=r,nu=nu);

u_rms = 0.05;
u_max = pi*u_rms;
fdFlux = SingleModeForcedDissipativeQGPVEMasked(wvt,k_f=k_f,k_r=k_r,u_rms=u_rms);
model = WaveVortexModel(wvt,nonlinearFlux=fdFlux);

% [u,v] = wvt.velocityField();
% U = max(max(max( sqrt(u.*u + v.*v) )))

deltaT = (wvt.x(2)-wvt.x(1))*0.25/u_max;
model.setupIntegrator(deltaT=deltaT, outputInterval=86400);

% model.integrateOneTimeStep();

% model.integrateToNextOutputTime();
model.integrateToTime(500*86400);

u2 = wvt.u.^2 + wvt.v.^2;
u_max = max(sqrt(u2(:)))
u_rms = sqrt(mean(u2(:)))

% figure, pcolor(wvt.x,wvt.y,wvt.ssh.'), shading interp

% model.integrateToNextOutputTime();
A02 = (wvt.A0_TE_factor/wvt.h) .* (wvt.A0.*conj(wvt.A0));
u_rms_alt = sqrt(2*sum(A02(:)))
GeostrophicEnergyK = wvt.transformToRadialWavenumber(A02);
dk = wvt.kRadial(2)-wvt.kRadial(1);
figure
plot(wvt.kRadial,GeostrophicEnergyK/dk), xlog, ylog
ylabel('m^3/s^2')
xlabel('1/m')
vlines(k_f,'g--')
% 
% figure, pcolor(log10(abs(fftshift(fdFlux.damp)))), shading interp
% figure, pcolor(log10(abs(fftshift(fdFlux.F)))), shading interp
% 
% DampK = wvt.transformToRadialWavenumber(fdFlux.damp);