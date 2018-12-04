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
k = 10.^linspace(log10(2*pi/1e6),log10(2*pi/10),150);
dk = diff(k);

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
xlim([0 1.1*max(1e4*Euv_const)])

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

Lx = 800e3;

figure('Name', 'Variances vs wavenumber')
subplot(3,1,1)
plot(k,k.*squeeze(Suv_const(50,:))), xlog, hold on, vlines(2*pi/Lx,'g--')
subplot(3,1,2)
plot(k,k.*squeeze(Seta_const(50,:))), xlog, hold on, vlines(2*pi/Lx,'g--')
subplot(3,1,3)
plot(k,k.*squeeze(Sw_const(50,:))), xlog, hold on, vlines(2*pi/Lx,'g--')

[omega2, k2, m2] = GMConst.SquaredFrequencyForWavenumber(k);

% summed by *mode* --- if summed by vertical wavenumber, need dm
[Suv_kj,j] = GMConst.HorizontalVelocitySpectrumAtWavenumberAndMode(k);
[Seta_kj,j] = GMConst.IsopycnalSpectrumAtWavenumberAndMode(k);
[Sw_kj,j] = GMConst.VerticalVelocitySpectrumAtWavenumberAndMode(k);

HorizontalShearNonlinearity = Suv_kj.*(k2./omega2);
IsopycnalSlopeNonlinearity = Seta_kj.*k2;
GradientFroudeNonlinearity = Suv_kj.*(m2./(N0*N0)); % Richardson number
VerticalStrainNonlinearity = Seta_kj.*m2;
return
m = j*pi/L;

% This is probably the right thing to do, because its the integral that
% sets the total energy of the wave.
VariancePreserving = (k').*j;

% Measures of stokes drift
StokesDriftF2 = Suv_kj.^2.*k2./omega2;
figure('Name', 'Stokes drift vs mode and wavenumber'), jpcolor(k,j,(VariancePreserving.*abs(StokesDriftF2)).'); shading flat, xlog, ylog
figure('Name', 'Stokes drift vs wavenumber'), plot(k,(k').*sum(StokesDriftF2,2)), xlog, hold on, vlines(2*pi/Lx,'g--')

figure('Name', 'Nonlinearity vs modes and wavenumber')
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

figure('Name', 'Variances vs modes and wavenumber')
subplot(3,1,1)
jpcolor(k,j,((k').*j.*abs(Suv_kj)).'); shading flat, xlog, ylog
subplot(3,1,2)
jpcolor(k,j,((k').*j.*abs(Seta_kj)).'); shading flat, xlog, ylog
subplot(3,1,3)
jpcolor(k,j,((k').*j.*abs(Sw_kj)).'); shading flat, xlog, ylog

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Stokes drift as function of wavenumber and depth
%
z = reshape(z,[],1);
k = reshape(k,1,[]);
m=shiftdim(j*pi/L,-1);
dk = cat(2,k(2)-k(1),diff(k));

% confirm that our summation actually works how we expected
U2 = dk.*shiftdim(Suv_kj,-1).*cos(z.*m).^2;
figure, plot(sum(sum(U2,3),2)*1e4,z);
xlabel('cm^2/s^2'), ylabel('depth')
title('horizontal velocity variance')

StokesDriftF2=dk.^2 .* shiftdim(Suv_kj.^2.*k2./omega2/4,-1);

Stokes = sum(sum(StokesDriftF2.*cos(2*z.*m).^2,3),2);
figure, plot(Stokes*1e4,z);
xlabel('cm^2/s^2'), ylabel('depth')
title('horizontal stokes drift variance')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Group velocity
%
% want k x m
cg_h = k.*m.^2*(N0^2-f0^2)./((k.^2 + m.^2).^(3/2) .* sqrt(N0^2*k.^2 + f0^2.*m.^2));

% figure('Name', 'Group velocity')
% jpcolor(k,j,cg_h.'); shading flat, xlog, ylog

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Decorrelation time should be 1./(k*cg_h)
%

T_decorrelation = 1./(cg_h .* k);
% T_decorrelation = max(T_decorrelation,[],3);

kappa = sum(sum(StokesDriftF2.*cos(2*z.*m).^2.*T_decorrelation,3),2);
figure, plot(kappa,z)
xlabel('m^2/s'), ylabel('depth')