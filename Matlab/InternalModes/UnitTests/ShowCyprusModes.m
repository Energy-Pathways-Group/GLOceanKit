load('CyprusSimulationParams.mat');

z = linspace(-Lz,0,nz).';
latitude = asind(coriolis/( 2*(7.2921e-5)));

im = InternalModes(rhobar,z,z,latitude,'method','finiteDifference'); % 'spectral'
im.upperBoundary = UpperBoundary.rigidLid;
im.lowerBoundary = LowerBoundary.freeSlip;
[F,G,h] = im.ModesAtWavenumber(0);
kappa = InternalModes.ConditionNumberAsFunctionOfModeNumber(F);

figure, plot(kappa), ylog

im = InternalModes(rhobar,z,z,latitude,'method','spectral'); % 'spectral'
im.upperBoundary = UpperBoundary.rigidLid;
im.lowerBoundary = LowerBoundary.freeSlip;
[F,G,h] = im.ModesAtWavenumber(0);
kappa = InternalModes.ConditionNumberAsFunctionOfModeNumber(F);

hold on
plot(kappa)

% 
% figure, plot(F(:,1:4),z)
% im.ShowLowestModesAtWavenumber(0);
