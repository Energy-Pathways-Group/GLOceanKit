if 1 == 0
    [rho, N2, zIn] = InternalModes.StratificationProfileWithName('constant');
    z = linspace(min(zIn),max(zIn),5000);
    N0 = sqrt(max(N2(z)));
else
    GM = GarrettMunkSpectrum('constant');
    z = GM.zInternal;
    N0 = GM.N_max;
end

if ~exist('GM','var')
    GM = GarrettMunkSpectrum(rho,[-L 0],latitude);
end

L_gm = 1300;
f0 = GM.f0;

Euv = GM.HorizontalVelocityVariance(z);
Eeta = GM.IsopycnalVariance(z);
Ew = GM.VerticalVelocityVariance(z);

N2 = GM.N2(z);
N = sqrt(N2);

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
plot(1e4*(Euv + Ew + N2.*Eeta)/2,z)
xlabel('total energy (cm^2/s^2)')
ylabel('depth (m)')

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
Siso = GM.IsopycnalSpectrumAtFrequencies([0 -L_gm/2 -L_gm],omega);
plot(omega,Siso), ylog, xlog
Sref = omega.^(-2); Sref(omega<f0) = 0; refIndex = find(omega>f0,1,'first'); Sref = Sref * (Siso(2,refIndex)/Sref(refIndex))*10;
hold on, plot(omega,Sref,'k','LineWidth',2)
ylim([1e1 3e6])
xlim(1.05*[0 N0])
title('horizontal isopycnal spectra')
xlabel('radians per second')

subplot(2,2,4)
omega = linspace(0,N0,500);
Sw = GM.VerticalVelocitySpectrumAtFrequencies(linspace(-500,0,20),omega);
plot(omega,Sw), ylog, xlog
ylim([1e-6 1e-1])
xlim(1.05*[0 N0])
title('horizontal vertical velocity spectra')
xlabel('radians per second')


figure
plot(GM.IsopycnalVariance(z),z)

GMConst = GarrettMunkSpectrumConstantStratification(N0,GM.z_in,GM.latitude);
[S, m] = GM.IsopycnalSpectrumAtVerticalWavenumbers();
[S_const, m_const] = GMConst.IsopycnalSpectrumAtVerticalWavenumbers();
figure
plot(m,S),ylog, xlog, hold on
plot(m_const,S_const)

trapz(z,GMConst.IsopycnalVariance(z))/GM.Lz
sum(S_const)*m_const(1)
sum(S)*m(1)

return

k = linspace(0,pi/10,150)';
S = GM.HorizontalVelocitySpectrumAtWavenumbers(k);

S( S<1e-2 ) = 1e-2;
figure, plot(k,S), ylog, xlog
ylim([1e-2 1.1*max(max(S))])
