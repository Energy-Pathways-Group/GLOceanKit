upperBoundary = UpperBoundary.freeSurface;
normalization = Normalization.kConstant;
nPoints = 100;

[rhoFunc, ~, zIn] = InternalModes.StratificationProfileWithName('exponential');
z = linspace(min(zIn),max(zIn),1024)';
im = InternalModesSpectral(rhoFunc,zIn,z,33,'nEVP',512);
im.normalization = normalization;
im.upperBoundary = upperBoundary;

z_g = im.GaussQuadraturePointsForModesAtFrequency(nPoints,0);

im_gauss = InternalModesSpectral(rhoFunc,zIn,z_g,33,'nEVP',512,'nModes',nPoints);
im_gauss.normalization = normalization;
im_gauss.upperBoundary = upperBoundary;

[F,G,h,~,uMaxRatio,wMaxRatio] = im_gauss.ModesAtFrequency(0,'uMax','wMax');
[F_hr,G_hr] = im.ModesAtFrequency(0);

if im_gauss.upperBoundary ==  UpperBoundary.rigidLid
    maxGMode = nPoints-2;
    figure('Name','Gauss-quadrature points, rigid lid')
elseif im_gauss.upperBoundary == UpperBoundary.freeSurface
    maxGMode = nPoints-1;
    figure('Name','Gauss-quadrature points, free surface')
end
subplot(1,2,1)
plot(G_hr(:,maxGMode)*wMaxRatio(maxGMode),im.z), hold on
scatter(G(:,maxGMode)*wMaxRatio(maxGMode),z_g)
title(sprintf('Highest resolvable G mode (%d)',maxGMode))
subplot(1,2,2)
plot(F_hr(:,maxGMode+1)*uMaxRatio(maxGMode+1),im.z), hold on
scatter(F(:,maxGMode+1)*uMaxRatio(maxGMode+1),z_g)
title(sprintf('Highest resolvable F mode (%d)',maxGMode+1))

kappaF = InternalModes.ConditionNumberAsFunctionOfModeNumberForModeIndices(F.*uMaxRatio,1:nPoints);
kappaG = InternalModes.ConditionNumberAsFunctionOfModeNumberForModeIndices(G.*wMaxRatio,1:nPoints);

return

fprintf('The condition number for the first %d modes of F is %f\n', maxGMode+1, cond(F(:,1:(maxGMode+1))));
fprintf('The condition number for the first %d modes of G is %f\n', maxGMode, cond(G(:,1:maxGMode)));

gamma = zeros(1,size(G,2));
for iMode = 1:size(G,2)
    gamma(iMode) = norm(G(:,iMode));
end
gamma = median(gamma)./gamma;
G_tilde = G .* gamma;
return
if im_gauss.upperBoundary == UpperBoundary.freeSurface
    gamma = norm(G(:,2))/norm(G(:,1));
    G = cat(2, gamma*G(:,1),G(:,2:end));
    F = cat(2, gamma*F(:,1),F(:,2:end));
    fprintf('After rescaling-------------\n')
    fprintf('The condition number for the first %d modes of F is %f\n', maxGMode+1, cond(F(:,1:(maxGMode+1))));
    fprintf('The condition number for the first %d modes of G is %f\n', maxGMode, cond(G(:,1:maxGMode)));
end

fprintf('Note that the barotropic mode destroys the conditioning. If you leave it off, everything does well again\n');

