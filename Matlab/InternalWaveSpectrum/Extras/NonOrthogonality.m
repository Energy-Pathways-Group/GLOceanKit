profilename = 'latmix-site1';
% profilename = 'exponential';
boundaryCondition = UpperBoundary.rigidLid;
latitude = 31;

[rho, N2, zIn] = InternalModes.StratificationProfileWithName(profilename);

z = linspace(min(zIn),max(zIn),5000)';

im = InternalModesAdaptiveSpectral(rho,zIn,z,latitude);

k = 2*pi/1000;

[F1,G1,h1,omega1] = im.ModesAtWavenumber(k);
[F2,G2,h2,omega2] = im.ModesAtWavenumber(2*k);

% mode 2 from k and mode 4 from 2k have similar frequencies
[F,G,h,k] = im.ModesAtFrequency(omega2(4));

figure
subplot(1,2,1)
plot(F1(:,2),z), hold on, plot(F2(:,4),z)
subplot(1,2,2)
plot(F(:,[2,4]),z)
