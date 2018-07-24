upperBoundary = UpperBoundary.rigidLid;
nPoints = 64;
k = 0;

[rhoFunc, ~, zIn] = InternalModes.StratificationProfileWithName('constant');
z = linspace(min(zIn),max(zIn),1024)';
im = InternalModesSpectral(rhoFunc,zIn,z,33,'nEVP',512);
im.upperBoundary = upperBoundary;

z_g = im.GaussQuadraturePointsForModesAtWavenumber(nPoints,k);

im_gauss = InternalModesSpectral(rhoFunc,zIn,z_g,33,'nEVP',256,'nModes',nPoints);
im_gauss.upperBoundary = upperBoundary;

[F,G,h] = im_gauss.ModesAtWavenumber(k);
[F_hr,G_hr] = im.ModesAtWavenumber(k);

figure('Name','Gauss-quadrature points')
if im_gauss.upperBoundary ==  UpperBoundary.rigidLid
    maxGMode = nPoints-2;
elseif im_gauss.upperBoundary == UpperBoundary.freeSurface
    maxGMode = nPoints-1;
end
subplot(1,2,1)
plot(G_hr(:,maxGMode),im.z), hold on
scatter(G(:,maxGMode),z_g)
title('Highest resolvable G mode')
subplot(1,2,2)
plot(F_hr(:,maxGMode+1),im.z), hold on
scatter(F(:,maxGMode+1),z_g)
title('Highest resolvable F mode')


fprintf('The condition number for the first %d modes of F is %f\n', maxGMode+1, cond(F(:,1:(maxGMode+1))));
fprintf('The condition number for the first %d modes of G is %f\n', maxGMode, cond(G(:,1:maxGMode)));

fprintf('Note that the barotropic mode destroys the conditioning. If you leave it off, everything does well again');