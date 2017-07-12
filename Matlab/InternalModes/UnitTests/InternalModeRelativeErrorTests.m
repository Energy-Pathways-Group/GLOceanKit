
methods = cell(5,1);
methods{1} = 'finiteDifference';
im = InternalModes('constant', 'finiteDifference', 64);


im.upperBoundary = 'free_surface';
im.normalization = 'const_F_norm';

im.ShowRelativeErrorAtFrequency(im.f0);
% im.ShowRelativeErrorAtWavenumber(0.1);