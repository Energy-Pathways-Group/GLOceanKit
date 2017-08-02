methods = cell(5,1);
methods{1} = 'densitySpectral';
methods{2} = 'wkbSpectral';
methods{3} = 'finiteDifference';
methods{4} = 'spectral';
im = InternalModes('constant', methods{1} , 64);

im.upperBoundary = 'free_surface';
im.normalization = 'const_G_norm';

N0 = 5.2e-3;
f0 = 7.9431e-05;
g = 9.81;
Lz = 5000;
k_star = sqrt( (N0*N0 - f0*f0) / (g*Lz) );

% im.ShowLowestModesAtWavenumber(0.1*k_star);
% return

im.ShowRelativeErrorAtWavenumber(0.1*k_star);
im.ShowRelativeErrorAtWavenumber(k_star);
im.ShowRelativeErrorAtWavenumber(10*k_star);

im.ShowRelativeErrorAtFrequency(0.1*N0);
im.ShowRelativeErrorAtFrequency(N0);
im.ShowRelativeErrorAtFrequency(10*N0);