profiles = cell(1,1);
profiles{1} = 'constant';
profiles{2} = 'exponential';

methods = cell(4,1);
methods{1} = 'finiteDifference';
methods{2} = 'wkbSpectral';
methods{3} = 'densitySpectral';
methods{4} = 'spectral';
methods{5} = 'wkbAdaptiveSpectral';

upperBoundary = UpperBoundary.freeSurface;

iProfile = 1;
iMethod = 1;

fprintf('**************************************************************\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the analytical solution
n = 2*64;
latitude = 33;
[rhoFunction, N2Function, zIn] = InternalModes.StratificationProfileWithName(profiles{iProfile});
z = linspace(min(zIn),max(zIn),n)';
if strcmp(profiles{iProfile},'constant')==1
    imAnalytical = InternalModesConstantStratification(5.2e-3,zIn,z,latitude,'nModes',n);
else
    imAnalytical = InternalModesExponentialStratification([5.2e-3 1300],zIn,z,latitude,'nModes',n);
end
imAnalytical.upperBoundary = upperBoundary;
imAnalytical.normalization = Normalization.kConstant;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute k_star, the wavenumber at which the free surface solution changes
N0 = 5.2e-3;
f0 = 7.9431e-05;
g = 9.81;
Lz = max(zIn)-min(zIn);
k_star = sqrt( (N0*N0 - f0*f0) / (g*Lz) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the error function: y is the true solution, x is the approximated
errorFunction = @(x,y) max(abs(x-y),[],1)./max(abs(y),[],1);
errorTolerance = 1e-2;

rhoInit = rhoFunction;
zInit = zIn;

im = InternalModes(rhoInit,zInit,z,latitude,'nModes',n, 'method', methods{iMethod});
im.upperBoundary = upperBoundary;
im.normalization = Normalization.kConstant;

k = 0.1*k_star;
[F,G,h,~,F2,N2G2,G2,uMaxRatio,wMaxRatio,kConstantRatio,omegaConstantRatio] = im.ModesAtWavenumber( k, 'F2', 'N2G2', 'G2', 'uMax', 'wMax', 'kConstant', 'omegaConstant' );
[F_analytical,G_analytical,h_analytical,~,F2_analytical,N2G2_analytical,G2_analytical,uMaxRatio_analytical,wMaxRatio_analytical,kConstantRatio_analytical,omegaConstantRatio_analytical] = imAnalytical.ModesAtWavenumber( k, 'F2', 'N2G2', 'G2', 'uMax', 'wMax', 'kConstant', 'omegaConstant' );

figure
subplot(2,2,1)
plot(F(:,1:4),im.z)
subplot(2,2,2)
plot(G(:,1:4),im.z)
subplot(2,2,3)
plot(F_analytical(:,1:4),im.z)
subplot(2,2,4)
plot(G_analytical(:,1:4),im.z)