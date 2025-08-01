%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Create a WVTransform
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lz = 4000;
N0 = 3*2*pi/3600; % buoyancy frequency at the surface, radians/seconds
L_gm = 1300; % thermocline exponential scale, meters
N2 = @(z) N0*N0*exp(2*z/L_gm);

wvt = WVTransformHydrostatic([800e3, 400e3, Lz], [256 128 40], N2=N2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Add a Garrett-Munk wave spectrum
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% wvt.initWithGMSpectrum(shouldShowDiagnostics=1);
wvt.initWithAlternativeSpectrum;


%%

Euv = squeeze(mean(mean(wvt.u.^2 + wvt.v.^2)));
Eeta = squeeze(mean(mean(wvt.eta.^2)));
Ew = squeeze(mean(mean(wvt.w.^2)));
E = 0.5*squeeze(mean(mean(wvt.u.^2 + wvt.v.^2 + shiftdim( wvt.N2,-2).*wvt.eta.^2)));

figure
tiledlayout(1,4,TileSpacing="tight")
nexttile
plot(1e4*Euv,wvt.z,LineWidth=2)
title('$E\left\langle u^2 + v^2 \right\rangle$','Interpreter','latex')
xlabel('cm^2/s^2')
ylabel('depth (m)')
xlim([0 1.05*max(1e4*Euv)])

nexttile
plot(Eeta,wvt.z,LineWidth=2)
title('$E\left\langle \eta^2 \right\rangle$','Interpreter','latex')
xlabel('m^2')
xlim([0 1.05*max(Eeta)])
set(gca,'YTickLabel',[]);

nexttile
plot(1e4*Ew,wvt.z,LineWidth=2)
title('$E\left\langle w^2 \right\rangle$','Interpreter','latex')
xlabel('cm^2/s^2')
xlim([0 1.05*max(1e4*Ew)])
set(gca,'YTickLabel',[]);

nexttile
plot(1e4*E,wvt.z,LineWidth=2)
title('$\frac{1}{2} E\left\langle u^2 + v^2 + N^2 \eta^2 \right\rangle$','Interpreter','latex')
xlabel('cm^2/s^2')
xlim([0 1.05*max(1e4*E)])
set(gca,'YTickLabel',[]);

fprintf('The total hydrostatic energy in the water column is %f J/m^2, compared to 3800 J/m^2 expected for GM.\n',wvt.rho0*trapz(wvt.z,E));
