%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CheckGMSpectrumFromModelOutput
%
% This scripts reads the model output generated from
% SaveModelOutputToNetCDF.m and does a few basics tests to see if the wave
% field is as expected.
%
% Jeffrey J. Early
% jeffrey@jeffreyearly.com
%
% December 6th, 2016      Version 1.0

file = '/Users/jearly/Desktop/InternalWaveModel_64x64x33@2017-03-28T120412.nc';
file = '/Users/jearly/Desktop/InternalWaveModel_512x64x33@2017-03-28T114337.nc';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Read in the problem dimensions and parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = ncread(file, 'x');
y = ncread(file, 'y');
z = ncread(file, 'z');
t = ncread(file, 't');

latitude = ncreadatt(file, '/', 'latitude');
N0 = ncreadatt(file, '/', 'N0');
f0 = 2*(7.2921e-5)*sin(latitude*pi/180);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Check the frequency spectrum at a given depth, but multiple (x,y)
% locations.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

depth = 1250;
[depth_index] = find(z <= depth, 1, 'last');

% This is what is wonderful about NetCDF. We can pull out a slice in time
% at a given depth, for a range of (strided!) x,y values.
stride = 4;
t_index = length(t)-1;
u3d = double(squeeze(ncread(file, 'u', [1 1 depth_index 1], [length(x)/stride length(y)/stride 1 t_index], [stride stride 1 1])));
v3d = double(squeeze(ncread(file, 'v', [1 1 depth_index 1], [length(x)/stride length(y)/stride 1 t_index], [stride stride 1 1])));

% Create a few 'mooring' time series from this.
[M, N, K] = size(u3d);
cv_mooring = zeros([K 1]);
subsample = 1;
iMooring = 0;
for i=1:subsample:M
	for j=1:subsample:N
		iMooring = iMooring+1;
		cv_mooring(:,iMooring) = squeeze(u3d(i,j,:) + sqrt(-1)*v3d(i,j,:));
	end
end

taper_bandwidth = 2;
psi=[];
%  [psi,lambda]=sleptap(size(cv_mooring,1),taper_bandwidth);
[omega_p, Spp, Snn, Spn] = mspec(t(2)-t(1),cv_mooring,psi);

omega = [ -flipud(omega_p(2:end)); omega_p];
[S_gm] = GarrettMunkHorizontalKineticEnergyRotarySpectrumWKB( omega, latitude, N0, 0 );
% We want the integral of this to give us the variance back, so we need to
% divide by 2*pi
S = (1/(2*pi))*[flipud(vmean(Snn,2)); vmean(Spp(2:end,:),2)];
figure, plot(omega,S), ylog
hold on, plot(omega,S_gm*0.226)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Check the variance of u,v,w,zeta as a function of depth
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uvVariance = zeros(length(z),1);
zetaVariance = zeros(length(z),1);
wVariance = zeros(length(z),1);
timeEnsemble = 1:6:min(400,length(t)-1); % apparently this is irrelevant. why?
for iTime = timeEnsemble
    u = double(squeeze(ncread(file, 'u', [1 1 1 iTime], [length(x) length(y) length(z) 1], [1 1 1 1])));
    v = double(squeeze(ncread(file, 'v', [1 1 1 iTime], [length(x) length(y) length(z) 1], [1 1 1 1])));
    w = double(squeeze(ncread(file, 'w', [1 1 1 iTime], [length(x) length(y) length(z) 1], [1 1 1 1])));
    zeta = double(squeeze(ncread(file, 'zeta', [1 1 1 iTime], [length(x) length(y) length(z) 1], [1 1 1 1])));

    uvVariance = uvVariance + squeeze(mean(mean(u.*u + v.*v,1),2));
    zetaVariance = zetaVariance + squeeze(mean(mean(zeta.*zeta,1),2));
    wVariance = wVariance + squeeze(mean(mean(w.*w,1),2));
end
uvVariance = uvVariance/length(timeEnsemble);
zetaVariance = zetaVariance/length(timeEnsemble);
wVariance = wVariance/length(timeEnsemble);


j_star = 3;
L_gm = 1.3e3; % thermocline exponential scale, meters
invT_gm = 5.2e-3; % reference buoyancy frequency, radians/seconds
E_gm = 6.3e-5; % non-dimensional energy parameter
E = L_gm*L_gm*L_gm*invT_gm*invT_gm*E_gm;
N2 = N0*N0;
B0 = 1/(pi/2 - atan( f0/sqrt(N2-f0*f0)));
H0 = 1/sum((j_star+(1:1024)).^(-5/2));

Phi = zeros(length(z),1);
Gamma = zeros(length(z),1);
D = max(abs(z)); g=9.81; N2=N0*N0; f0 = 2*(7.2921e-5)*sin(latitude*pi/180);
for j=1:1024
   Phi = Phi + (2/D)*H0*((j_star+j).^(-5/2)) * cos(z*j*pi/D).^2;
   Gamma = Gamma + (2/D)*H0*((j_star+j).^(-5/2)) * sin(z*j*pi/D).^2;
end

ExpectedUV = E*Phi*( 3*N2/2 - f0*f0 - (B0/2)*f0*sqrt(N2-f0*f0) )/(N2-f0*f0);
ExpectedZeta = E*Gamma*(1/2 - (B0/2)*(f0/N2)*sqrt(N2-f0*f0))/(N2-f0*f0);
ExpectedW = E*Gamma*f0*f0*((B0/f0)*sqrt(N2-f0*f0)-1)/(N2-f0*f0);

ExpectedTotal = ExpectedUV + ExpectedW + N2*ExpectedZeta;

totalVariance = uvVariance + wVariance + N0*N0.*zetaVariance;

figure
subplot(1,4,1)
% plot([44 44], [z(1) z(end)], 'k' ,'LineWidth', 2)
plot(1e4*uvVariance,z,'LineWidth', 2), hold on
plot(1e4*ExpectedUV,z,'LineWidth', 2)
xlabel('cm^2/s^2'), ylabel('depth (m)')
title(sprintf('HKE (%.2fGM)',trapz(z,uvVariance)/trapz(z,ExpectedUV) ))

subplot(1,4,2)
% plot([44 44], [z(1) z(end)], 'k' ,'LineWidth', 2)
plot(1e4*wVariance,z,'LineWidth', 2), hold on
plot(1e4*ExpectedW,z,'LineWidth', 2)
xlabel('cm^2/s^2'), set(gca,'YTickLabel',[]);
title(sprintf('VKE (%.2fGM)',trapz(z,wVariance)/trapz(z,ExpectedW) ))

subplot(1,4,3)
% plot([53 53], [z(1) z(end)], 'k' ,'LineWidth', 2), hold on
plot(zetaVariance,z,'LineWidth', 2), hold on
plot(ExpectedZeta,z,'LineWidth', 2)
xlabel('m^2'), set(gca,'YTickLabel',[]);
title(sprintf('isopycnal variance (%.2fGM)',trapz(z,zetaVariance)/trapz(z,ExpectedZeta) ))

subplot(1,4,4)
% plot([30 30], [z(1) z(end)], 'k' ,'LineWidth', 2), hold on
plot(1e4*totalVariance ,z,'LineWidth', 2), hold on
plot(1e4*ExpectedTotal ,z,'LineWidth', 2)
xlabel('cm^2 s^{-2}'), set(gca,'YTickLabel',[]);
title(sprintf('total energy (%.2fGM)',trapz(z,totalVariance)/trapz(z,ExpectedTotal) ))
