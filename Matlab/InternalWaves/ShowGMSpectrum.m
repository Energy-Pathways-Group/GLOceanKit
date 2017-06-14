latitude = 33;
f0 = 2 * 7.2921E-5 * sin( latitude*pi/180 );
rho0 = 1025;
g = 9.81;
N0 = 5.2e-3/4; % 'average' GM bouyancy frequency, sqrt( (1-exp(-L/L_gm))/(exp(L/L_gm-1) -1))
rho = @(z) -(N0*N0*rho0/g)*z + rho0;

N0 = 5.2e-3;
L_gm = 1.3e3;
rho = @(z) rho0*(1 + L_gm*N0*N0/(2*g)*(1 - exp(2*z/L_gm)));
L = 5000;


z = linspace(-L,0,100)';

% S = GarrettMunkHorizontalVelocitySpectrum( omega, latitude, rho, [-L 0], zOut,'exact' );
% S( S<1e-2 ) = 1e-2;
% figure, plot(omega,S), ylog
% ylim([1e0 1.1*max(max(S))])

GM = GarrettMunkSpectrum(rho,[-L 0], z,latitude);

% omega = linspace(-N0,N0,200);
[E,S,omega] = GM.HorizontalVelocityVariance();
figure, plot(1e4*E,z)
xlabel('velocity variance (cm^2/s^2)')
ylabel('depth (m)')

return

k = linspace(0,pi/10,150)';
S = GM.HorizontalVelocitySpectrumAtWavenumbers(k);

S( S<1e-2 ) = 1e-2;
figure, plot(k,S), ylog, xlog
ylim([1e-2 1.1*max(max(S))])
