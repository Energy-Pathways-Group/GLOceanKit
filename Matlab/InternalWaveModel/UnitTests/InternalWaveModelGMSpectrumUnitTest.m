%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% InternalWaveModelGMSpectrumUnitTest
%
% This script uses the InternalWaveModel to create, and validate, a
% Garrett-Munk spectrum in a linear internal wave field with constant
% stratification.
%
% Jeffrey J. Early
% jeffrey@jeffreyearly.com
%
% March 25th, 2016      Version 1.0
% March 30th, 2016      Version 1.1
% November 17th, 2016   Version 1.2
% February 22nd, 2017   Version 1.3


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lx = 800e3;
Ly = 100e3;
Lz = 5000;

Nx = 512;
Ny = 512/8;
Nz = 65; % 2^n + 1 grid points, to match the Winters model, but 2^n ok too.

latitude = 31;
N0 = 5.2e-3; % Choose your stratification 7.6001e-04

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wavemodel = InternalWaveModel([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);
wavemodel.InitializeWithGMSpectrum(1.0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Time step the model and save a time series of (u,v) from the center of
% the domain.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

depth = Lz/2;
depthIndex = find(wavemodel.z-Lz > -depth,1,'first');
stride = 4;
t = 0:15*60:4*86400;
xIndices = 1:stride:Nx;
yIndices = 1:stride:Ny;
cv_mooring = zeros([length(t) length(xIndices)*length(yIndices)]);
startTime = datetime('now');
for iTime=1:length(t)
    if mod(iTime,10) == 0
        timePerStep = (datetime('now')-startTime)/(iTime-1);
        timeRemaining = (length(t)-iTime+1)*timePerStep;   
        fprintf('\ttime step %d of %d. Estimated finish time %s (%s from now)\n', iTime, length(t), datestr(datetime('now')+timeRemaining), datestr(timeRemaining, 'HH:MM:SS')) ;
    end
    
    [u,v]=wavemodel.VelocityFieldAtTime(t(iTime));
    
    cv_mooring(iTime,:) = reshape(u(xIndices,yIndices,depthIndex),1,[]) + sqrt(-1)*reshape(v(xIndices,yIndices,depthIndex),1,[]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot the rotary spectrum of the horizontal velocity
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[omega_p, Spp, Snn, Spn] = mspec(t(2)-t(1),cv_mooring,[]);
omega = [ -flipud(omega_p(2:end)); omega_p];

% We want the integral of this to give us the variance back, so we need to
% divide by 2*pi
S = (1/(2*pi))*[flipud(vmean(Snn,2)); vmean(Spp(2:end,:),2)];
figure, semilogy(omega*86400/(2*pi),S, 'LineWidth', 2)

% Compare it to the non-hydrostatic version of the GM spectrum
L_GM = 1.3e3; % thermocline exponential scale, meters
N_GM = 5.2e-3; % reference buoyancy frequency, radians/seconds
E_GM = 6.3e-5; % non-dimensional energy parameter
f0 = wavemodel.f0;
nonHydrostaticFactor = ((N0*N0 - omega.*omega)/(N0*N0));
prefactor = abs(real(L_GM*L_GM*N0*N_GM*E_GM/pi)*(f0./omega).*(1./sqrt(omega.*omega-f0*f0).*nonHydrostaticFactor));
S_gm = real(prefactor.*((f0 - omega).^2)./(omega.*omega));
S_gm(abs(omega)<f0) = 0;

hold on, plot(omega*86400/(2*pi),S_gm, 'k' ,'LineWidth', 2)
xlabel('frequency (cycles per day)'), ylabel('m^2/s^2'), title('horizontal velocity power spectrum')
legend('model output', 'GM')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot the structure of the vertical variances
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = 0;
[u,v]=wavemodel.VelocityFieldAtTime(t);
[w,zeta] = wavemodel.VerticalFieldsAtTime(t);

z = wavemodel.z;
uvVariance = squeeze(mean(mean(u.*u + v.*v,1),2));
zetaVariance = squeeze(mean(mean(zeta.*zeta,1),2));

zeta2 = squeeze(mean(mean(zeta.*zeta,1),2));
u2 = squeeze(mean(mean(u.*u,1),2)+mean(vmean(v.*v,1),2));
w2 = squeeze(mean(mean(w.*w,1),2));

figure
subplot(1,3,1)
plot([44 44], [z(1) z(end)], 'k' ,'LineWidth', 2), hold on
plot(1e4*uvVariance,z,'LineWidth', 2)
xlabel('cm^2/s^2'), ylabel('depth (m)'), title('horizontal velocity variance')
xlim([0 1.05*max(1e4*uvVariance)])
legend('GM', 'model output')

subplot(1,3,2)
plot([53 53], [z(1) z(end)], 'k' ,'LineWidth', 2), hold on
plot(zetaVariance,z,'LineWidth', 2)
xlabel('m^2'), set(gca,'YTickLabel',[]), title('isopycnal variance')

subplot(1,3,3)
plot([30 30], [z(1) z(end)], 'k' ,'LineWidth', 2), hold on
plot(0.5*1e4*(u2 + w2 + N0*N0.*zeta2) ,z,'LineWidth', 2)
xlabel('cm^2 s^{-2}'), set(gca,'YTickLabel',[]); title('total energy'),
xlim([0 1.05*max(0.5*1e4*(u2 + w2 + N0*N0.*zeta2))])

HKE = 0.5*(u.*u + v.*v);
HKE_int = trapz(z,HKE,3);

VKE = 0.5*(w.*w);
VKE_int = trapz(z,VKE,3);

PE = 0.5*(N0^2)*zeta.*zeta;
PE_int = trapz(z,PE,3);

totalGM = mean(mean(HKE_int + VKE_int + PE_int))*1032; % scaled by the density of water
fprintf('The total energy in the water column is %f J/m^2, compared to 3800 J/m^2 expected for GM. This is %.2f times larger.\n',totalGM,totalGM/3800);
fprintf('We expect the integrated energy to be D/b time larger than GM, or %.1f\n',Lz/1300);

w = round(Nz/4);
c = round(Nz/2);
indices=(c-w):(c+w);
AvgHKE = mean(mean(mean(HKE(:,:,indices))))*1e4;
fprintf('The average 2*HKE is %.2f cm^2/s^2 in the middle half of the domain, compared to 44 cm^2/s^2 expected for GM.\n',2*AvgHKE);