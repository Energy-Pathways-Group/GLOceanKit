n = 64;
lat = 33;
N0 = 5.2e-3; % reference buoyancy frequency, radians/seconds
g = 9.81;
rho_0 = 1025;
zIn = [-5000 0];
zOut = linspace(zIn(1),0,2000)';
L_gm = 1.3e3; % thermocline exponential scale, meters
rho = @(z) rho_0*(1 + L_gm*N0*N0/(2*g)*(1 - exp(2*z/L_gm)));
im = InternalModesAdaptiveSpectral(rho,zIn,zOut,lat);

omega = 3*im.f0;
[F,G,h] = im.ModesAtFrequency( 0.95*N0 );

figure
subplot(1,3,1)
plot(F(:,1:4),im.z, 'LineWidth', 2)
ylabel('depth (meters)');
xlabel('(u,v)-modes');

b = subplot(1,3,2);
plot(G(:,1:4),im.z, 'LineWidth', 2)
xlabel('w-modes');
ytick([]);

subplot(1,3,3)
plot(sqrt(im.N2),im.z, 'LineWidth', 2), hold on
% if ~isempty(im.N2Function)
%     plot(sqrt(im.N2Function(self.z)),im.z, 'LineWidth', 2)
% end
xlim([0.0 1.1*max(sqrt(im.N2))])
xlabel('buoyancy frequency');
ytick([]);