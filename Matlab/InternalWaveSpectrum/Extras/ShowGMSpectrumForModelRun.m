latitude = 33;
f0 = 2 * 7.2921E-5 * sin( latitude*pi/180 );
rho0 = 1025;
g = 9.81;

N0 = 5.2e-3;
L_gm = 1.3e3;
rho = @(z) rho0*(1 + L_gm*N0*N0/(2*g)*(1 - exp(2*z/L_gm)));
L = 4000;

z = linspace(-L,0,100)';

if ~exist('GM','var')
    GM = GarrettMunkSpectrum(rho,[-L 0],latitude);
end
Euv = GM.HorizontalVelocityVariance(z);
Eeta = GM.IsopycnalVariance(z);
Ew = GM.VerticalVelocityVariance(z);
N2 = GM.N2(z);
N = sqrt(N2);

L_const = L/2;
z_const = linspace(-L_const,0,100)';

N0_const = N0/sqrt(2);
GMConst = GarrettMunkSpectrumConstantStratification(N0_const,[-L_const 0],latitude);
EnergyScale = L/L_gm/2.5;
Euv_const = EnergyScale*GMConst.HorizontalVelocityVariance(z_const);
Eeta_const = EnergyScale*GMConst.IsopycnalVariance(z_const);
Ew_const = EnergyScale*GMConst.VerticalVelocityVariance(z_const);

figure
subplot(1,3,1)
plot(1e4*Euv,z), hold on
plot(1e4*Euv_const,z_const)
xlabel('velocity variance (cm^2/s^2)')
ylabel('depth (m)')
subplot(1,3,2)
plot(Eeta,z), hold on
plot(Eeta_const,z_const)
xlabel('isopycnal variance (m^2)')
ylabel('depth (m)')
subplot(1,3,3)
plot(1e4*(Euv + Ew + N2.*Eeta)/2,z), hold on
plot(1e4*(Euv_const + Ew_const + N0_const*N0_const.*Eeta_const)/2,z_const)
xlabel('total energy (cm^2/s^2)')
ylabel('depth (m)')
legend('Exponential profile (GM reference)', 'Constant stratification')

figure
subplot(1,2,1)
plot(1e4*Euv.*(N0./N),z), hold on
plot(1e4*Euv_const.*(N0./N0_const),z_const)
vlines(44,'k--')
xlabel('velocity variance (cm^2/s^2)')
ylabel('depth (m)')
title('wkb scaled')
subplot(1,2,2)
plot(Eeta.*(N/N0),z),hold on
plot(Eeta_const*(N0_const/N0),z_const)
vlines(53,'k--')
xlabel('isopycnal variance (m^2)')
ylabel('depth (m)')
title('wkb scaled')

E = 0.5*(Euv + Ew + N2.*Eeta);
Etotal = trapz(z,E);

figure

subplot(2,1,1)
omega = linspace(-N0,N0,200);
S = GM.HorizontalVelocitySpectrumAtFrequencies([0 -L_gm/2 -L_gm],omega);
S_const = GMConst.HorizontalVelocitySpectrumAtFrequencies([0 -L_gm/2 -L_gm],omega);;
plot(omega,S), ylog, hold on,
ax = gca;
ax.ColorOrderIndex = 1;
plot(omega,S_const,'--')
ylim([1e-4 1e2])
xlim(1.05*[-N0 N0])
title('horizontal velocity spectra')
xlabel('radians per second')
legend('surface','650 m', '1300 m')

subplot(2,1,2)
omega = linspace(0,N0,500);
Siso = GM.IsopycnalSpectrumAtFrequencies([0 -L_gm/2 -L_gm],omega);
Siso_const = GMConst.IsopycnalSpectrumAtFrequencies([0 -L_gm/2 -L_gm],omega);
plot(omega,Siso), ylog, xlog, hold on,
ax = gca;
ax.ColorOrderIndex = 1;
plot(omega,Siso_const,'--')
Sref = omega.^(-2); Sref(omega<f0) = 0; refIndex = find(omega>f0,1,'first'); Sref = Sref * (Siso(2,refIndex)/Sref(refIndex))*10;
hold on, plot(omega,Sref,'k','LineWidth',2)
ylim([1e1 3e6])
xlim(1.05*[0 N0])
title('isopycnal spectra')
xlabel('radians per second')
legend('surface','650 m', '1300 m')

% subplot(2,2,4)
% omega = linspace(0,N0,500);
% Sw = GM.HorizontalVerticalVelocitySpectrumAtFrequencies(linspace(-500,0,20),omega);
% plot(omega,Sw), ylog, xlog
% ylim([1e-6 1e-1])
% xlim(1.05*[0 N0])
% title('horizontal vertical velocity spectra')
% xlabel('radians per second')


return

k = linspace(0,pi/10,150)';
S = GM.HorizontalVelocitySpectrumAtWavenumbers(k);

S( S<1e-2 ) = 1e-2;
figure, plot(k,S), ylog, xlog
ylim([1e-2 1.1*max(max(S))])
