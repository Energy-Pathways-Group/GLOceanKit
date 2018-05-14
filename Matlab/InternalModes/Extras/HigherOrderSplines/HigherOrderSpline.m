[rhoFunc, ~, zIn] = InternalModes.StratificationProfileWithName('exponential');
z = linspace(min(zIn),max(zIn),64)';

K=4;
z_knot = NaturalKnotsForSpline( z, K );
B = bspline( z, z_knot, K );
X = squeeze(B(:,:,1));
% N = size(X,1);
% M = size(X,2);
rho = rhoFunc(z);
% want rho = X*m
m = X\rho;

% im = InternalModesSpectral(rhoFunc,zIn,z,33,'nEVP',128);

