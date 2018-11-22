latitude = 33;
f0 = 2 * 7.2921E-5 * sin( latitude*pi/180 );
rho0 = 1025;
g = 9.81;
N0 = 5.2e-3; % 'average' GM bouyancy frequency, sqrt( (1-exp(-L/L_gm))/(exp(L/L_gm-1) -1))
rho = @(z) -(N0*N0*rho0/g)*z + rho0;

L_gm = 1.3e3;
L = 1300;

z = linspace(-L,0,100)';

GMConst = GarrettMunkSpectrumConstantStratification(N0,[-L 0],latitude);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The primary variances as a function of depth
%
Euv_const = GMConst.HorizontalVelocityVariance(z);
Eeta_const = GMConst.IsopycnalVariance(z);
Ew_const = GMConst.VerticalVelocityVariance(z);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Same variances, but now computed by summing their associated spectra
%
k = 10.^linspace(log10(2*pi/1e6),log10(2*pi/10),150); dk = k(2)-k(1);

Suv_const = GMConst.HorizontalVelocitySpectrumAtWavenumbers(z,k);
Seta_const = GMConst.IsopycnalSpectrumAtWavenumbers(z,k);
Sw_const = GMConst.VerticalVelocitySpectrumAtWavenumbers(z,k);

Euv_const_summed = trapz(k,Suv_const,2);
Eeta_const_summed = trapz(k,Seta_const,2);
Ew_const_summed = trapz(k,Sw_const,2);

figure
subplot(1,3,1)
plot(1e4*Euv_const,z,'r'),hold on
plot(1e4*Euv_const_summed,z,'r--')
legend('analytical','analytical summed')
xlabel('horizontal velocity variance (cm^2/s^2)')
ylabel('depth (m)')
xlim([0 1.1*max(1e4*Euv)])

subplot(1,3,2)
plot(Eeta_const,z,'r'),hold on
plot(Eeta_const_summed,z,'r--')
legend('analytical','analytical summed')
xlabel('isopycnal variance (m^2)')
ylabel('depth (m)')

subplot(1,3,3)
plot(1e4*Ew_const,z,'r'),hold on
plot(1e4*Ew_const_summed,z,'r--')
legend('analytical','analytical summed')
xlabel('vertical velocity variance (cm^2/s^2)')
ylabel('depth (m)')

Lx = 108e3;

figure
subplot(3,1,1)
plot(k,k.*squeeze(Suv_const(50,:))), xlog, hold on, vlines(2*pi/Lx,'g--')
subplot(3,1,2)
plot(k,k.*squeeze(Seta_const(50,:))), xlog, hold on, vlines(2*pi/Lx,'g--')
subplot(3,1,3)
plot(k,k.*squeeze(Sw_const(50,:))), xlog, hold on, vlines(2*pi/Lx,'g--')

[omega2, k2, m2] = GMConst.SquaredFrequencyForWavenumber(k);

[Suv_kj,j] = GMConst.HorizontalVelocitySpectrumAtWavenumberAndMode(k);
[Seta_kj,j] = GMConst.IsopycnalSpectrumAtWavenumberAndMode(k);
[Sw_kj,j] = GMConst.VerticalVelocitySpectrumAtWavenumberAndMode(k);

HorizontalShearNonlinearity = Suv_kj.*(k2./omega2);
IsopycnalSlopeNonlinearity = Seta_kj.*k2;
GradientFroudeNonlinearity = Suv_kj.*(m2./(N0*N0)); % Richardson number
VerticalStrainNonlinearity = Seta_kj.*m2;

% This is probably the right thing to do, because its the integral that
% sets the total energy of the wave.
VariancePreserving = (k').*j;

% Measures of stokes drift
StokesDriftF2 = Suv_kj.^2.*k2./omega2;
figure, jpcolor(k,j,(VariancePreserving.*abs(StokesDriftF2)).'); shading flat, xlog, ylog
figure, plot(k,(k').*sum(StokesDriftF2,2)), xlog

figure
subplot(4,1,1)
jpcolor(k,j,(VariancePreserving.*abs(HorizontalShearNonlinearity)).'); shading flat, xlog, ylog
colorbar
subplot(4,1,2)
jpcolor(k,j,(VariancePreserving.*abs(IsopycnalSlopeNonlinearity)).'); shading flat, xlog, ylog
colorbar
subplot(4,1,3)
jpcolor(k,j,(VariancePreserving.*abs(GradientFroudeNonlinearity)).'); shading flat, xlog, ylog
colorbar
subplot(4,1,4)
jpcolor(k,j,(VariancePreserving.*abs(VerticalStrainNonlinearity)).'); shading flat, xlog, ylog
colorbar

figure
subplot(3,1,1)
jpcolor(k,j,((k').*j.*abs(Suv_kj)).'); shading flat, xlog, ylog
subplot(3,1,2)
jpcolor(k,j,((k').*j.*abs(Seta_kj)).'); shading flat, xlog, ylog
subplot(3,1,3)
jpcolor(k,j,((k').*j.*abs(Sw_kj)).'); shading flat, xlog, ylog