latitude = 33;

[rho, ~, zIn] = InternalModes.StratificationProfileWithName('constant');
L = max(zIn)-min(zIn);

if ~exist('GM','var')
    GM = GarrettMunkSpectrum(rho,[-L 0],latitude);
end
GMConst = GarrettMunkSpectrumConstantStratification(5.2e-3,zIn,latitude);

k = linspace(0,pi/10,150)';
S = GM.HorizontalVelocitySpectrumAtWavenumbers(-2500,k);
S = GM.HorizontalVelocitySpectrumAtWavenumbers(-2500,GM.k);
figure, plot(GM.k,S), ylog, xlog

S = GMConst.HorizontalVelocitySpectrumAtWavenumbers(-2500,GM.k);
hold on, plot(GM.k,S)

figure,plot(GM.k(2:end),diff(log10(S))./diff(log10(GM.k))), xlog

return

z = GM.zInternal;
k = GM.k;
j=shiftdim(1:GM.nModes,-1);

N2 = GM.N_max*GM.N_max;
f2 = GM.f0*GM.f0;
f = GM.f0;

m = j*pi/GM.Lz;
omega2 = (N2-f2)*(k.^2./(k.^2 + m.^2)) + f2;

% Phi analytical and numerical are spot on.
Phi = (2/GM.Lz)*GM.H(j) .* (m.^2./(k.^2 + m.^2)) .* cos(z.*m).^2;
Bfunc = (2/pi)*((f*m.^2)./(N2*k.^2 + f2*m.^2)) .* sqrt( (N2-f2)./(k.^2 + m.^2));
C = squeeze(1+f2./omega2);

Nmax = GM.N_max;
C_function = @(omega) (abs(omega)<f | abs(omega) > Nmax)*0 + (abs(omega) >= f & abs(omega) <= Nmax).*( (1+(f./omega).^2) );
C2 = C_function(GM.omega_k);

return

Euv = GM.HorizontalVelocityVariance(z);
Eeta = GM.HorizontalIsopycnalVariance(z);
Ew = GM.HorizontalVerticalVelocityVariance(z);
N2 = GM.N2(z);
N = sqrt(N2);

N0_const = N0/2;
GMConst = GarrettMunkSpectrumConstantStratification(N0_const,[-L 0],latitude);
EnergyScale = L/L_gm/2;
Euv_const = EnergyScale*GMConst.HorizontalVelocityVariance(z);
Eeta_const = EnergyScale*GMConst.HorizontalIsopycnalVariance(z);
Ew_const = EnergyScale*GMConst.HorizontalVerticalVelocityVariance(z);

figure
subplot(1,3,1)
plot(1e4*Euv,z), hold on
plot(1e4*Euv_const,z)
xlabel('velocity variance (cm^2/s^2)')
ylabel('depth (m)')
subplot(1,3,2)
plot(Eeta,z), hold on
plot(Eeta_const,z)
xlabel('isopycnal variance (m^2)')
ylabel('depth (m)')
subplot(1,3,3)
plot(1e4*(Euv + Ew + N2.*Eeta)/2,z), hold on
plot(1e4*(Euv_const + Ew_const + N0_const*N0_const.*Eeta_const)/2,z)
xlabel('total energy (cm^2/s^2)')
ylabel('depth (m)')
legend('Exponential profile (GM reference)', 'Constant stratification')

figure
subplot(1,2,1)
plot(1e4*Euv.*(N0./N),z), hold on
plot(1e4*Euv_const.*(N0./N0_const),z)
vlines(44,'k--')
xlabel('velocity variance (cm^2/s^2)')
ylabel('depth (m)')
title('wkb scaled')
subplot(1,2,2)
plot(Eeta.*(N/N0),z),hold on
plot(Eeta_const*(N0_const/N0),z)
vlines(53,'k--')
xlabel('isopycnal variance (m^2)')
ylabel('depth (m)')
title('wkb scaled')

E = 0.5*(Euv + Ew + N2.*Eeta);
Etotal = trapz(z,E);

return

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
