[rhoFunc, ~, zIn] = InternalModes.StratificationProfileWithName('exponential');
z = linspace(min(zIn),max(zIn),1024)';
im = InternalModesSpectral(rhoFunc,zIn,z,33,'nEVP',512);

z_g = im.GaussQuadraturePointsForModesAtWavenumber(64,0);

im_gauss = InternalModesSpectral(rhoFunc,zIn,z_g,33,'nEVP',256,'nModes',64);
[F,G,h] = im_gauss.ModesAtWavenumber(0);

cond(F)
cond((im_gauss.N2-im_gauss.f0^2).*G)