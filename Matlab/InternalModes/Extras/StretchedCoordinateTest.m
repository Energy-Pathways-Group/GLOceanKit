[rhoFunction, N2Function, zIn] = InternalModes.StratificationProfileWithName('exponential');
n = 64;
latitude = 33;
z = linspace(min(zIn),max(zIn),n)';
im = InternalModes(rhoFunction,zIn,z,latitude,'nModes',n, 'method', 'wkbSpectral');


% x = chebfun(N2Function,zIn);
% y = inv(cumsum(sqrt(x)));