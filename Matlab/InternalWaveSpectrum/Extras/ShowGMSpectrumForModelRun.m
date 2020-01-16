shouldSavePlots = 0;

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

N2 = GM.N2(z);
N = sqrt(N2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize a stratification with good parameters for the model
%
% L_const = 1300;
% z_const = linspace(-L_const,0,100)';
% N0_const = N0;
% GMConst = GarrettMunkSpectrumConstantStratification(N0_const,[-L_const 0],latitude);
% EnergyScale = 1;

omega = 2*(2*pi/86400);
N0_const = N0/4;
L_const = 1300*sqrt( (N0*N0 + omega*omega)/(N0_const*N0_const - omega*omega));
z_const = linspace(-L_const,0,100)';
GMConst = GarrettMunkSpectrumConstantStratification(N0_const,[-L_const 0],latitude);
EnergyScale = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Look at the primary variances versus depth
%
Euv = GM.HorizontalVelocityVariance(z);
Eeta = GM.IsopycnalVariance(z);
Ew = GM.VerticalVelocityVariance(z);

Euv_const = EnergyScale*GMConst.HorizontalVelocityVariance(z_const);
Eeta_const = EnergyScale*GMConst.IsopycnalVariance(z_const);
Ew_const = EnergyScale*GMConst.VerticalVelocityVariance(z_const);

ExpLineStyle = '-';
ConstLineStyle = '--';

figure('Name','Variances Vs Depth')
subplot(1,4,1)
plot(1e4*Euv,z,ExpLineStyle,'LineWidth',2,'Color','k'), hold on
plot(1e4*Euv_const,z_const,ConstLineStyle,'LineWidth',2,'Color','k')
xlabel('cm^2/s^2')
ylabel('depth (m)')
title('$\left\langle u^2 \right\rangle$','interpreter','latex')

subplot(1,4,2)
plot(1e4*Ew,z,ExpLineStyle,'LineWidth',2,'Color','k'), hold on
plot(1e4*Ew_const,z_const,ConstLineStyle,'LineWidth',2,'Color','k')
xlabel('cm^2/s^2')
title('$\left\langle w^2 \right\rangle$','interpreter','latex')

subplot(1,4,3)
plot(Eeta,z,ExpLineStyle,'LineWidth',2,'Color','k'), hold on
plot(Eeta_const,z_const,ConstLineStyle,'LineWidth',2,'Color','k')
xlabel('m^2')
title('$\left\langle \eta^2 \right\rangle$','interpreter','latex')

subplot(1,4,4)
plot(1e4*(Euv + Ew + N2.*Eeta)/2,z,ExpLineStyle,'LineWidth',2,'Color','k'), hold on
plot(1e4*(Euv_const + Ew_const + N0_const*N0_const.*Eeta_const)/2,z_const,ConstLineStyle,'LineWidth',2,'Color','k')
xlabel('cm^2/s^2')
title('$\left\langle u^2 + w^2 + N^2 \eta^2 \right\rangle$','interpreter','latex')

legend('Exponential', 'Constant','Location','southeast')
packfig(1,4)

if shouldSavePlots == 1
    print('-depsc2','VariancesVsDepth.eps')
end


figure('Name','WKBScaled Variances Vs Depth')
subplot(1,3,1)
plot(1e4*Euv.*(N0./N),z,ExpLineStyle,'LineWidth',2,'Color','k'), hold on
plot(1e4*Euv_const.*(N0./N0_const),z_const,ConstLineStyle,'LineWidth',2,'Color','k')
plot(44*ones(size(z)),z,'LineWidth',2,'Color',0.5*[1 1 1])
xlabel('cm^2/s^2')
ylabel('depth (m)')
title('$\frac{N_0}{N(z)}\left\langle u^2 \right\rangle$','interpreter','latex')

subplot(1,3,2)
plot(1e4*Ew.*(N/N0),z,ExpLineStyle,'LineWidth',2,'Color','k'),hold on
plot(1e4*Ew_const*(N0_const/N0),z_const,ConstLineStyle,'LineWidth',2,'Color','k')
xlabel('cm^2/s^2')
ylabel('depth (m)')
title('$\frac{N(z)}{N_0}\left\langle w^2 \right\rangle$','interpreter','latex')

subplot(1,3,3)
plot(Eeta.*(N/N0),z,ExpLineStyle,'LineWidth',2,'Color','k'),hold on
plot(Eeta_const*(N0_const/N0),z_const,ConstLineStyle,'LineWidth',2,'Color','k')
plot(53*ones(size(z)),z,'LineWidth',2,'Color',0.5*[1 1 1])
xlabel('m^2')
ylabel('depth (m)')
title('$\frac{N(z)}{N_0}\left\langle \eta^2 \right\rangle$','interpreter','latex')

packfig(1,3)

legend('Exponential', 'Constant','GM Ref','Location','southeast')

if shouldSavePlots == 1
    print('-depsc2','VariancesVsDepth_WKBScaled.eps')
end

E = 0.5*(Euv + Ew + N2.*Eeta);
Etotal = trapz(z,E);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Look at frequency spectra of the same quantities
%

omega_hke = linspace(-N0,N0,200);
omega_iso = linspace(f0,N0,4000);
depths = [-100 -325 -650];
legendlabels = cell(length(depths),1);
for iDepth = 1:length(depths)
   legendlabels{iDepth} = sprintf('%d m',abs(depths(iDepth)));
end

Suv = GM.HorizontalVelocitySpectrumAtFrequencies(depths,omega_hke);
Siso = GM.IsopycnalSpectrumAtFrequencies(depths,omega_iso);
Sw = GM.VerticalVelocitySpectrumAtFrequencies(depths,omega_iso);

Suv_const = EnergyScale*GMConst.HorizontalVelocitySpectrumAtFrequencies(depths,omega_hke);
Siso_const = EnergyScale*GMConst.IsopycnalSpectrumAtFrequencies(depths,omega_iso);
Sw_const = EnergyScale*GMConst.VerticalVelocitySpectrumAtFrequencies(depths,omega_iso);

frequencyscale = 3600/(2*pi);

figure('Name','HKE Frequency Spectrum')
plot(frequencyscale*omega_hke,Suv,ExpLineStyle,'LineWidth',2), ylog, hold on,
ax = gca;
ax.ColorOrderIndex = 1;
plot(frequencyscale*omega_hke,Suv_const,ConstLineStyle,'LineWidth',2)
ylim([1e-4 20])
xlim(1.05*frequencyscale*[-N0 N0])
title('$\left\langle u^2 \right\rangle$(f)','interpreter','latex')
xlabel('cph')
legend(legendlabels)

if shouldSavePlots == 1
    print('-depsc2','HorizontalVelocityFrequencySpectrum.eps')
end

figure('Name','Spectra')
subplot(2,1,1)
plot(frequencyscale*omega_iso,Siso,ExpLineStyle,'LineWidth',2), ylog, xlog, hold on,
ax = gca;
ax.ColorOrderIndex = 1;
plot(frequencyscale*omega_iso,Siso_const,ConstLineStyle,'LineWidth',2)
% Sref = omega_iso.^(-2); Sref(omega_iso<f0) = 0; refIndex = find(omega_iso>f0,1,'first'); Sref = Sref * (Siso(2,refIndex)/Sref(refIndex))*10;
% hold on, plot(omega_iso,Sref,'k','LineWidth',2)
ylim([3e2 3e6])
xlim(frequencyscale*[f0 1.05*N0])
title('$\left\langle \eta^2 \right\rangle$(f)','interpreter','latex')
xlabel('cph')
legend(legendlabels)

subplot(2,1,2)
plot(frequencyscale*omega_iso,Sw,ExpLineStyle,'LineWidth',2), ylog, xlog, hold on,
ax = gca;
ax.ColorOrderIndex = 1;
plot(frequencyscale*omega_iso,Sw_const,ConstLineStyle,'LineWidth',2)
ylim([2e-4 1e-1])
xlim(frequencyscale*[f0 1.05*N0])
title('$\left\langle w^2\right\rangle$(f) ','interpreter','latex')
xlabel('cph')

if shouldSavePlots == 1
    print('-depsc2','IsopycnalAndVerticalVelocityFrequencySpectra.eps')
end

packfig(2,1)

return

% subplot(2,2,4)
% omega = linspace(0,N0,500);
% Sw = GM.HorizontalVerticalVelocitySpectrumAtFrequencies(linspace(-500,0,20),omega);
% plot(omega,Sw), ylog, xlog
% ylim([1e-6 1e-1])
% xlim(1.05*[0 N0])
% title('horizontal vertical velocity spectra')
% xlabel('radians per second')

[S,m] = GM.IsopycnalSpectrumAtVerticalWavenumbers();
[S_const,m_const] = GMConst.IsopycnalSpectrumAtVerticalWavenumbers();

figure
plot(m,S),  xlog, ylog, hold on
plot(m_const,S_const)

return

k = linspace(0,pi/10,150)';
S = GM.HorizontalVelocitySpectrumAtWavenumbers(k);

S( S<1e-2 ) = 1e-2;
figure, plot(k,S), ylog, xlog
ylim([1e-2 1.1*max(max(S))])
