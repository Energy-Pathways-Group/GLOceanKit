% We are ultimately going to use one grid, so it may as well be the
% quadrature grid of the hydrostatic modes.
% 
% I think we need to compute the inverse transformation in its quadrature
% grid, then project it back onto physical space. So, if Tq is (zq x m) and
% if V (zq x m) is the eigenmatrix in physical space, then T^{-1} V is the
% matrix in Chebyshev space. Or T \hat{V}, inv(T \hat{V}) =
% inv(\hat{V})*T^{-1}
% 
% Thus, we need to determine
% 1) 
[rhoFunc, ~, zIn] = InternalModes.StratificationProfileWithName('exponential');
zIn = [-4000, 0];
z = linspace(min(zIn),max(zIn),1024)';
upperBoundary = UpperBoundary.rigidLid;
normalization = Normalization.kConstant;
nPoints = 10;

% Find the quadrature points for k=0
k = 0;
im = InternalModesSpectral(rho=rhoFunc,zIn=zIn,zOut=z,latitude=33,nEVP=256);
im.normalization = normalization;
im.upperBoundary = upperBoundary;
z = im.GaussQuadraturePointsForModesAtWavenumber(nPoints,k);

% Compute F,G at those points
im = InternalModesSpectral(rho=rhoFunc,zIn=zIn,zOut=z,latitude=33,nEVP=256);
im.normalization = normalization;
im.upperBoundary = upperBoundary;
[Finv,Ginv,self.h] = im.ModesAtWavenumber(k);



wavelength = 1e3;
k = (2*pi)./wavelength;


for iK=1:length(k)
    

    im = InternalModesSpectral(rho=rhoFunc,zIn=zIn,zOut=z,latitude=33,nEVP=256);
    im.normalization = normalization;
    im.upperBoundary = upperBoundary;
    
    z_g(:,iK) = im.GaussQuadraturePointsForModesAtWavenumber(nPoints,k(iK));
end