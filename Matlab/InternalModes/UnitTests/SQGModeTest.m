n = 1000;

lat = 33;
N0 = 5.2e-3; % reference buoyancy frequency, radians/seconds
g = 9.81;
rho_0 = 1025;
zIn = [-5000 0];
z = linspace(zIn(1),0,n)';

rho = @(z) -(N0*N0*rho_0/g)*z + rho_0;
N2 = @(z) N0*N0*ones(size(z));

im = InternalModesConstantStratification( rho, zIn, z, lat );

k = 10.^linspace(log10(1e-5),log10(1e-1),10)';
psi_t = im.SurfaceModesAtWavenumber( k );
psi_b = im.BottomModesAtWavenumber( k );

im_spec = InternalModesFiniteDifference( rho, zIn, z, lat );
psi_t_spec = im_spec.SurfaceModesAtWavenumber( k );
psi_b_spec = im_spec.BottomModesAtWavenumber( k );

im_wkb = InternalModesWKB( rho, zIn, z, lat );
psi_t_wkb = im_wkb.SurfaceModesAtWavenumber(k(1));

figure
subplot(2,1,1)
plot(im.f0*squeeze(psi_t),z)
hold on, plot(im.f0*psi_t_spec,z)
subplot(2,1,2)
plot(im.f0*squeeze(psi_b),z)
hold on, plot(im.f0*psi_b_spec,z)

[rho, N2Func, zIn] = InternalModes.StratificationProfileWithName('exponential');
im_spec = InternalModesSpectral( rho, zIn, z, lat );
im_wkb = InternalModesWKB( rho, zIn, z, lat );

psi_t_spec = im_spec.SurfaceModesAtWavenumber( k(1) );
psi_t_wkb = im_wkb.SurfaceModesAtWavenumber(k(1));
figure
plot(im.f0*squeeze(psi_t_spec),z),
hold on
plot(im.f0*squeeze(psi_t_wkb),z)