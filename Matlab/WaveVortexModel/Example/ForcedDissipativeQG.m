%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lx = 50e3;
Ly = 50e3;

Nx = 256;
Ny = 256;

latitude = 25;

wvt = WaveVortexTransformSingleMode([Lx, Ly], [Nx, Ny], h=0.8, latitude=latitude);

% Forcing parameters
dk = (wvt.k(2)-wvt.k(1));
k_f = 30*dk;
k_r = 4*dk;
u_rms = 0.05;


fdFlux = SingleModeForcedDissipativeQGPVEMasked(wvt,k_f=k_f,k_r=k_r,u_rms=u_rms);
model = WaveVortexModel(wvt,nonlinearFlux=fdFlux);

u_max = pi*u_rms; % not sure how much bigger u_max will be than u_rms
deltaT = (wvt.x(2)-wvt.x(1))*0.25/u_max;
model.setupIntegrator(deltaT=deltaT, outputInterval=86400);

model.createNetCDFFileForModelOutput('ForcedDissipativeQG.nc',shouldOverwriteExisting=1);
model.integrateToTime(250*86400);

% model.integrateToNextOutputTime();

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
vlines([k_f,k_r],'g--')
% 
% figure, pcolor(log10(abs(fftshift(fdFlux.damp)))), shading interp
% figure, pcolor(log10(abs(fftshift(fdFlux.F)))), shading interp
% 
% DampK = wvt.transformToRadialWavenumber(fdFlux.damp);