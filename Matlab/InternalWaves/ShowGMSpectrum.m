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

GM = GarrettMunkSpectrum(rho,[-L 0],latitude);

Euv = GM.HorizontalVelocityVariance(z);
Eeta = GM.HorizontalIsopycnalVariance(z);
Ew = GM.HorizontalVerticalVelocityVariance(z);

figure
subplot(1,3,1)
plot(1e4*Euv,z)
xlabel('velocity variance (cm^2/s^2)')
ylabel('depth (m)')
subplot(1,3,2)
plot(Eeta,z)
xlabel('isopycnal variance (m^2)')
ylabel('depth (m)')
subplot(1,3,3)
plot(1e4*Ew,z)
xlabel('vertical velocity variance (cm^2/s^2)')
ylabel('depth (m)')

N2 = GM.N2(z);
N = sqrt(N2);
figure
subplot(1,2,1)
plot(1e4*Euv.*(N0./N),z), hold on
vlines(44,'k--')
xlabel('velocity variance (cm^2/s^2)')
ylabel('depth (m)')
title('wkb scaled')
subplot(1,2,2)
plot(Eeta.*(N/N0),z)
vlines(53,'k--')
xlabel('isopycnal variance (m^2)')
ylabel('depth (m)')
title('wkb scaled')

E = 0.5*(Euv + Ew + N2.*Eeta);
Etotal = trapz(z,E);

figure

subplot(2,2,[1 2])
omega = linspace(-N0,N0,200);
S = GM.HorizontalVelocitySpectrumAtFrequencies([0 -L_gm/2 -L_gm],omega);
plot(omega,S), ylog
ylim([1e-4 1e2])
xlim(1.05*[-N0 N0])
title('horizontal velocity spectra')
xlabel('radians per second')

subplot(2,2,3)
omega = linspace(0,N0,500);
Siso = GM.HorizontalIsopycnalSpectrumAtFrequencies([0 -L_gm/2 -L_gm],omega);
plot(omega,Siso), ylog, xlog
Sref = omega.^(-2); Sref(omega<f0) = 0; refIndex = find(omega>f0,1,'first'); Sref = Sref * (Siso(2,refIndex)/Sref(refIndex))*10;
hold on, plot(omega,Sref,'k','LineWidth',2)
ylim([1e1 3e6])
xlim(1.05*[0 N0])
title('horizontal isopycnal spectra')
xlabel('radians per second')

subplot(2,2,4)
omega = linspace(0,N0,500);
Sw = GM.HorizontalVerticalVelocitySpectrumAtFrequencies(linspace(-500,0,20),omega);
plot(omega,Sw), ylog, xlog
ylim([1e-6 1e-1])
xlim(1.05*[0 N0])
title('horizontal vertical velocity spectra')
xlabel('radians per second')


return

k = linspace(0,pi/10,150)';
S = GM.HorizontalVelocitySpectrumAtWavenumbers(k);

S( S<1e-2 ) = 1e-2;
figure, plot(k,S), ylog, xlog
ylim([1e-2 1.1*max(max(S))])
