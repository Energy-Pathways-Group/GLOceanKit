methods = cell(5,1);
methods{1} = 'finiteDifference';
methods{2} = 'wkbSpectral';
methods{3} = 'densitySpectral';
methods{4} = 'spectral';
methods{5} = 'wkbAdaptiveSpectral';

upperBoundary = UpperBoundary.rigidLid;
N0 = 5.2e-3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the analytical solution
n = 2*64;
latitude = 33;
[rhoFunction, N2Function, zIn] = InternalModes.StratificationProfileWithName('exponential');
z = linspace(min(zIn),max(zIn),n)';
imAnalytical = InternalModesExponentialStratification([N0 1300], zIn, z, latitude,'nModes',n);
imAnalytical.upperBoundary = upperBoundary;
imAnalytical.normalization = Normalization.kConstant;

% Define the error function: y is the true solution, x is the approximated
errorFunction = @(x,y) max(abs(x-y),[],1)./max(abs(y),[],1);
errorTolerance = 1e-2;

iMethod = 2;
if iMethod > 1
    im = InternalModes(rhoFunction,zIn,z,latitude,'nModes',n, 'method', methods{iMethod}, 'nEVP', n);
else
    im = InternalModes(rhoFunction,zIn,z,latitude,'nModes',n, 'method', methods{iMethod});
end
im.upperBoundary = upperBoundary;
im.normalization = Normalization.kConstant;

omega = 0.1*N0;
[F,G,h] = im.ModesAtFrequency( omega );
[F_analytical,G_analytical,h_analytical] = imAnalytical.ModesAtFrequency( omega );
max_error = max([errorFunction(h,h_analytical); errorFunction(F,F_analytical); errorFunction(G,G_analytical)],[],1);

iModes = 1;
figure
subplot(1,2,1)
plot(F(:,iModes),z), hold on,
plot(F_analytical(:,iModes),z)
subplot(1,2,2)
plot(G(:,iModes),z), hold on,
plot(G_analytical(:,iModes),z)