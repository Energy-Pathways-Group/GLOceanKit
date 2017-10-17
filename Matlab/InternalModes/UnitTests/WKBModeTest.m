%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Stratification
%
[rhoFunction, N2Function, zIn] = InternalModes.StratificationProfileWithName('exponential');
N0 = 5.2e-3;
zOut = linspace(min(zIn),max(zIn),500)';
latitude = 33;

im = InternalModesWKB(rhoFunction, zIn, zOut, latitude);
omega = im.f0 + 0.4*(N0-im.f0);
[F, G, h] = im.ModesAtFrequency(omega);


imAnalytical = InternalModesWKBSpectral(rhoFunction, zIn, zOut, latitude,'nEVP',512);
[F_analytical, G_analytical, h_analytical] = imAnalytical.ModesAtFrequency(omega);

iMode = 4;
figure
subplot(1,2,1)
plot(G(:,iMode),im.z), hold on
plot(G_analytical(:,iMode),im.z)
subplot(1,2,2)
plot(F(:,iMode),im.z), hold on
plot(F_analytical(:,iMode),im.z)