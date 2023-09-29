% The part of this that bothers me is the location of the bottom quadrature point
% for high wavenumbers. It sure seems like it should grab the extrema of
% the F mode at the bottom, but it doesn't. OTOH, we know that the F modes
% aren't the right thing.

[rhoFunc, ~, zIn] = InternalModes.StratificationProfileWithName('exponential');
zIn = [-4000, 0];
z = linspace(min(zIn),max(zIn),1024)';
upperBoundary = UpperBoundary.rigidLid;
normalization = Normalization.kConstant;
nPoints = 10;

wavelength = 10.^(5:-0.5:1);
k = (2*pi)./wavelength;
z_g = zeros(nPoints,length(k));
for iK=1:length(k)
    

    im = InternalModesSpectral(rho=rhoFunc,zIn=zIn,zOut=z,latitude=33,nEVP=512);
    im.normalization = normalization;
    im.upperBoundary = upperBoundary;
    
    z_g(:,iK) = im.GaussQuadraturePointsForModesAtWavenumber(nPoints,k(iK));
end

figure
subplot(1,3,1)
plot(sqrt(im.N2)*3600/(2*pi),im.z)
ylabel('depth')
xlabel('cph')
title('N(z)')

subplot(1,3,[2 3])
plot(log10(k),z_g)
xlabel('wavelength')
yticklabels([])
xlim([min(log10(k)) max(log10(k))])
xticks(log10(2*pi./[1e5 1e4 1e3 1e2 1e1]))
xticklabels({'100 km','10 km', '1 km', '100 m', '10 m'})
title('Quadrature points, mode 10')

% figure, plot(log10(wavelength),z_g)
% xlabel('wavelength (log10)')
% ylabel('depth')

return

im_gauss = InternalModesSpectral(rhoFunc,zIn,z_g,33,'nEVP',256,'nModes',nPoints);
im_gauss.normalization = normalization;
im_gauss.upperBoundary = upperBoundary;

[F,G,h] = im_gauss.ModesAtWavenumber(k);
[F_hr,G_hr] = im.ModesAtWavenumber(k);

if im_gauss.upperBoundary ==  UpperBoundary.rigidLid
    maxGMode = nPoints-2;
    figure('Name','Gauss-quadrature points, rigid lid')
elseif im_gauss.upperBoundary == UpperBoundary.freeSurface
    maxGMode = nPoints-1;
    figure('Name','Gauss-quadrature points, free surface')
end
subplot(1,2,1)
plot(G_hr(:,maxGMode),im.z), hold on
scatter(G(:,maxGMode),z_g)
title(sprintf('Highest resolvable G mode (%d)',maxGMode))
subplot(1,2,2)
plot(F_hr(:,maxGMode+1),im.z), hold on
scatter(F(:,maxGMode+1),z_g)
title(sprintf('Highest resolvable F mode (%d)',maxGMode+1))


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

