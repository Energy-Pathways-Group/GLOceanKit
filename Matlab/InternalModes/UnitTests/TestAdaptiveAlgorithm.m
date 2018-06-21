n = 64;
latitude = 33;
f0 = 2*(7.2921e-5)*sin(latitude*pi/180);

stratification = 'pycnocline-constant'; N0 = 0.0498; omega = f0 + 0.80*(N0-f0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Stratification
%
[rhoFunction, N2Function, zIn] = InternalModes.StratificationProfileWithName(stratification);

Lz=min(zIn);
z = linspace(Lz,0,n)';
rho = rhoFunction;
zOut = linspace(min(zIn),max(zIn),5000)';

im = InternalModesAdaptiveSpectral(rho,zIn,zOut,latitude);

[F,G,h] = im.ModesAtFrequency( omega );

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