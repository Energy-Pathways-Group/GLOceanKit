
% Fetch a pre-built profile from the internal modes class
[rho, N2, zIn] = InternalModes.StratificationProfileWithName('exponential');
latitude = 33;

% Initialize the GarrettMunkSpectrum class with the profile, but do *not*
% reinitialize it if it already exists. Recomputing the modes can take a
% while.
if ~exist('GM','var')
    GM = GarrettMunkSpectrum(rho,zIn,latitude);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute and plot the variances as a function of depth.
%
z = linspace(min(zIn),max(zIn),5000)';

approximation = 'gm';
approximation = 'exact';
Euv = GM.HorizontalVelocityVariance(z,approximation);
Eeta = GM.IsopycnalVariance(z,approximation);
Ew = GM.VerticalVelocityVariance(z,approximation);
E = 0.5*(Euv + Ew + N2(z).*Eeta); % total energy

fprintf('Total depth integrated energy: %f\n',trapz(z,E));

figure
subplot(1,4,1)
plot(1e4*Euv,z)
xlabel('velocity variance (cm^2/s^2)')
ylabel('depth (m)')
subplot(1,4,2)
plot(Eeta,z)
xlabel('isopycnal variance (m^2)')
set(gca,'YTickLabel',[]);
subplot(1,4,3)
plot(1e4*Ew,z)
xlabel('vertical velocity variance (cm^2/s^2)')
set(gca,'YTickLabel',[]);
subplot(1,4,4)
plot(1e4*E,z)
xlabel('total energy (cm^2/s^2)')
set(gca,'YTickLabel',[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the WKB scaled versions of the variances.
%

N0 = sqrt(max(N2(z)));
N = sqrt(GM.N2(z));

figure
subplot(1,3,1)
plot(1e4*Euv.*(N0./N),z), hold on
vlines(44,'k--')
xlabel('velocity variance (cm^2/s^2)')
ylabel('depth (m)')
title('wkb scaled')
legend('actual', 'gm ref','Location', 'northeast')
xlim([0 max(1e4*Euv.*(N0./N))*1.1]);

subplot(1,3,2)
plot(Eeta.*(N/N0),z)
vlines(53,'k--')
xlabel('isopycnal variance (m^2)')
set(gca,'YTickLabel',[]);
title('wkb scaled')
xlim([0 max(Eeta.*(N/N0))*1.1]);

subplot(1,3,3)
plot(1e4*Ew.*(N/N0),z)
wGM = GM.E*GM.f0*2/pi/GM.L_gm/GM.invT_gm;
vlines(1e4*wGM,'k--')
xlabel('vertical velocity variance (cm^2/s^2)')
set(gca,'YTickLabel',[]);
title('wkb scaled')
xlim([0 max(1e4*Ew.*(N/N0))*1.1]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the frequency spectra of the same quantities
%

L_gm = 1300;
f0 = GM.f0;

omega = linspace(-N0,N0,200);
depths = [0 -L_gm/2 -L_gm];
depthlabels = {sprintf('%d m',depths(1)), sprintf('%d m',depths(2)), sprintf('%d m',depths(3))};

Suv = GM.HorizontalVelocitySpectrumAtFrequencies(depths,omega,approximation);
Siso = GM.IsopycnalSpectrumAtFrequencies(depths,omega,approximation);
Sw = GM.VerticalVelocitySpectrumAtFrequencies(depths,omega,approximation);

figure

subplot(2,2,[1 2])
plot(omega,Suv), ylog
ylim([1e-4 1e2])
xlim(1.05*[-N0 N0])
title('horizontal velocity spectra')
xlabel('radians per second')
legend(depthlabels,'Location', 'northeast')

subplot(2,2,3)
plot(omega,Siso), ylog, xlog
Sref = omega.^(-2); Sref(omega<f0) = 0; refIndex = find(omega>f0,1,'first'); Sref = Sref * (Siso(2,refIndex)/Sref(refIndex))*10;
hold on, plot(omega,Sref,'k','LineWidth',2)
ylim([1e1 3e6])
xlim(1.05*[0 N0])
title('isopycnal spectra')
xlabel('radians per second')

subplot(2,2,4)

plot(omega,Sw), ylog, xlog
ylim([1e-6 1e-1])
xlim(1.05*[0 N0])
title('vertical velocity spectra')
xlabel('radians per second')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Work in progress
%
return
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
