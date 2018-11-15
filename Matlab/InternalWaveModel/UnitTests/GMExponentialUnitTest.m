%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GMExponentialUnitTest
%
% This script uses the InternalWaveModel to create, and validate, a
% Garrett-Munk spectrum in a linear internal wave field with exponential
% stratification.
%
% Jeffrey J. Early
% jeffrey@jeffreyearly.com
%
% November 5th, 2018      Version 1.0



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aspectRatio = 4;

Lx = 100e3;
Ly = aspectRatio*100e3;
Lz = 4000;

N = 32;
Nx = N;
Ny = aspectRatio*N;
Nz = N+1;

latitude = 31;
N0 = 5.2e-3;
b = 1300;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
% Initialization requires computing the internal modes, which can take
% quite a while. Therefore, we check to see if the model was previously
% initialized with the same parameters, and use that when possible.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('wavemodel','var') || any([wavemodel.Lx wavemodel.Ly wavemodel.Lz wavemodel.Nx wavemodel.Ny wavemodel.Nz wavemodel.N0 wavemodel.b] ~= [Lx Ly Lz Nx Ny Nz N0 b])
    wavemodel = InternalWaveModelExponentialStratification([Lx, Ly, Lz], [Nx, Ny, Nz], [N0 b], linspace(-Lz,0,Nz), latitude);
end
wavemodel.FillOutWaveSpectrum();
wavemodel.InitializeWithGMSpectrum(1.0);

[u,v,w,zeta]=wavemodel.VariableFieldsAtTime(0,'u','v','w','zeta');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute the spatially averaged variances from the dynamical variables
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uvVariance = squeeze(mean(mean(u.*u + v.*v,1),2));
zetaVariance = squeeze(mean(mean(zeta.*zeta,1),2));
wVariance = squeeze(mean(mean(w.*w,1),2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Let's compare this to the Garrett Munk spectrum
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('GM','var')
    GM = GarrettMunkSpectrum(wavemodel.internalModes.rhoFunction,[-Lz 0],latitude);
end
z = wavemodel.z;
Euv = GM.HorizontalVelocityVariance(z);
Eeta = GM.IsopycnalVariance(z);
Ew = GM.VerticalVelocityVariance(z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot everything
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(1,3,1)
plot([44 44], [z(1) z(end)], 'k' ,'LineWidth', 2), hold on
plot(1e4*(N0./sqrt(wavemodel.N2)).*uvVariance,z,'LineWidth', 2)
plot(1e4*(N0./sqrt(wavemodel.N2)).*GM.HorizontalVelocityVariance(z),z,'LineWidth', 2)
xlabel('cm^2/s^2'), ylabel('depth (m)'), title('horizontal velocity variance')
xlim([0 1.05*max(1e4*uvVariance)])
legend('GM', 'model output', 'theory')

subplot(1,3,2)
plot([53 53], [z(1) z(end)], 'k' ,'LineWidth', 2), hold on
plot((sqrt(wavemodel.N2)/N0).*zetaVariance,z,'LineWidth', 2)
plot((sqrt(wavemodel.N2)/N0).*GM.IsopycnalVariance(z),z,'LineWidth', 2)
xlabel('m^2'), set(gca,'YTickLabel',[]), title('isopycnal variance')

subplot(1,3,3)
plot([30 30], [z(1) z(end)], 'k' ,'LineWidth', 2), hold on
plot(0.5*1e4*(N0./sqrt(wavemodel.N2)).*(uvVariance + wVariance + wavemodel.N2.*zetaVariance) ,z,'LineWidth', 2)
xlabel('cm^2 s^{-2}'), set(gca,'YTickLabel',[]); title('total energy'),
xlim([0 1.05*max(0.5*1e4*(uvVariance + wVariance + wavemodel.N2.*zetaVariance))])
