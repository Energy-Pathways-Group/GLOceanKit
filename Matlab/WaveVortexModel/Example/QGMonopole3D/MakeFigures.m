wvt = WVTransform.waveVortexTransformFromFile('qg-eddy.nc',iTime=inf);

figure
tl = tiledlayout(3,1);
nexttile

ssh = wvt.seaSurfaceHeight;
pcolor(wvt.x/1e3, wvt.y/1e3, ssh.'), shading interp
axis equal

[~,I] = max(ssh(:));
[i1,i2] = ind2sub(size(ssh),I);


rho = wvt.rho_prime;
rho_total = wvt.rho_total;
sliceIndex = find(wvt.y< wvt.y(i2),1,'last');

nexttile
pcolor(wvt.x/1000,wvt.z,squeeze(rho(:,sliceIndex,:)).'); colorbar; clim([min(rho(:)),max(rho(:))]), shading interp, hold on
contour(wvt.x/1000,wvt.z,squeeze(rho_total(:,sliceIndex,:)).',linspace(min(rho_total(:)),max(rho_total(:)),10),'k','LineWidth',0.5);
xlabel('x (km)'), ylabel('z (m)')

nexttile
rv = wvt.diffX(wvt.v) - wvt.diffY(wvt.u);
pcolor(wvt.x/1000,wvt.z,squeeze(rv(:,sliceIndex,:)).'); colorbar; clim([-1,1]*max(abs(rv(:)))), shading interp
xlabel('x (km)'), ylabel('z (m)')