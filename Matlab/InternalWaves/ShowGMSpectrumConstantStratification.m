latitude = 33;
f0 = 2 * 7.2921E-5 * sin( latitude*pi/180 );
rho0 = 1025;
g = 9.81;
N0 = 5.2e-3; % 'average' GM bouyancy frequency, sqrt( (1-exp(-L/L_gm))/(exp(L/L_gm-1) -1))
rho = @(z) -(N0*N0*rho0/g)*z + rho0;

L_gm = 1.3e3;
L = 5000;

z = linspace(-L,0,100)';


GMConst = GarrettMunkSpectrumConstantStratification(N0,[-L 0],latitude);
Euv_const = GMConst.HorizontalVelocityVariance(z);
Eeta_const = GMConst.HorizontalIsopycnalVariance(z);
Ew_const = GMConst.HorizontalVerticalVelocityVariance(z);

omega = linspace(-N0,N0,2000); dOmega = omega(2)-omega(1);
Euv_const_summed = sum(GMConst.HorizontalVelocitySpectrumAtFrequencies(z,omega),2)*dOmega;
omega = linspace(0,N0,2000); dOmega = omega(2)-omega(1);
Eeta_const_summed = sum(GMConst.HorizontalIsopycnalSpectrumAtFrequencies(z,omega),2)*dOmega;
Ew_const_summed = sum(GMConst.HorizontalVerticalVelocitySpectrumAtFrequencies(z,omega),2)*dOmega;

GM = GarrettMunkSpectrum(rho,[-L 0],latitude);
Euv = GM.HorizontalVelocityVariance(z);
Eeta = GM.HorizontalIsopycnalVariance(z);
Ew = GM.HorizontalVerticalVelocityVariance(z);

N2 = GMConst.N2(z);
N = sqrt(N2);

figure
subplot(1,3,1)
plot(1e4*Euv,z), hold on
plot(1e4*Euv_const,z,'--')
plot(1e4*Euv_const_summed,z,'.')
legend('numerical','analytical','analytical summed')
xlabel('velocity variance (cm^2/s^2)')
ylabel('depth (m)')
subplot(1,3,2)
plot(Eeta,z), hold on
plot(Eeta_const,z,'--')
plot(Eeta_const_summed,z,'.')
legend('numerical','analytical','analytical summed')
xlabel('isopycnal variance (m^2)')
ylabel('depth (m)')
subplot(1,3,3)
plot(1e4*Ew,z), hold on
plot(1e4*Ew_const,z,'--')
plot(1e4*Ew_const_summed,z,'.')
legend('numerical','analytical','analytical summed')
xlabel('vertical velocity variance (cm^2/s^2)')
ylabel('depth (m)')

E = 0.5*(Euv + Ew + N2.*Eeta);
E_const = 0.5*(Euv_const + Ew_const + N2.*Eeta_const);
Etotal = trapz(z,E)
Etotal_const = trapz(z,E_const)

figure

subplot(2,2,[1 2])
omega = linspace(-N0,N0,2000);
Suv = GM.HorizontalVelocitySpectrumAtFrequencies([0 -L_gm/2 -L_gm],omega);
Suv_const = GMConst.HorizontalVelocitySpectrumAtFrequencies([0 -L_gm/2 -L_gm],omega);
plot(omega,Suv), ylog, hold on
plot(omega,Suv_const,'--'), ylog
ylim([1e-4 1e2])
xlim(1.05*[-N0 N0])
title('horizontal velocity spectra')
xlabel('radians per second')

subplot(2,2,3)
omega = linspace(0,N0,2000); 
Siso = GM.HorizontalIsopycnalSpectrumAtFrequencies([0 -L_gm/2 -L_gm],omega);
Siso_const = GMConst.HorizontalIsopycnalSpectrumAtFrequencies([0 -L_gm/2 -L_gm],omega);
plot(omega,Siso), ylog, xlog, hold on
plot(omega,Siso_const,'--'), ylog
ylim([1e1 3e6])
xlim(1.05*[0 N0])
title('horizontal isopycnal spectra')
xlabel('radians per second')

subplot(2,2,4)
Sw = GM.HorizontalVerticalVelocitySpectrumAtFrequencies([0 -L_gm/2 -L_gm],omega);
Sw_const = GMConst.HorizontalVerticalVelocitySpectrumAtFrequencies([0 -L_gm/2 -L_gm],omega);
plot(omega,Sw), ylog, xlog, hold on
plot(omega,Sw_const,'--'), ylog
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
