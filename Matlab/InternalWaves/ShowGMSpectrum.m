latitude = 33;
f0 = 2 * 7.2921E-5 * sin( latitude*pi/180 );
rho0 = 1025;
g = 9.81;
N0 = 5.2e-3/4; % 'average' GM bouyancy frequency, sqrt( (1-exp(-L/L_gm))/(exp(L/L_gm-1) -1))
rho = @(z) -(N0*N0*rho0/g)*z + rho0;
L = 5000;

omega = linspace(-2*N0,2*N0,100);
zOut = [-2500; -1250; 0];

S = GarrettMunkHorizontalVelocitySpectrum( omega, latitude, rho, [-L 0], zOut );

