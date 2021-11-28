% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults

file = '/Volumes/MoreStorage/Data/cyprus_eddy_wvm/cyprus_eddy-more-stratification-strong.nc';

% load('NonlinearFluxes.mat');
energies.t = ncread(file,'t');
energies.EnergyIGWPlus = ncread(file,'EnergyIGWPlus');
energies.EnergyIGWMinus = ncread(file,'EnergyIGWMinus');
energies.EnergyIOBaroclinic = ncread(file,'EnergyIOBaroclinic');
energies.EnergyIOBarotropic = ncread(file,'EnergyIOBarotropic');
energies.EnergyGeostrophicBaroclinic = ncread(file,'EnergyGeostrophicBaroclinic');
energies.EnergyGeostrophicBarotropic = ncread(file,'EnergyGeostrophicBarotropic');
% energies.EnergyResidual = ncread(file,'EnergyResidual');

energies.EnergyIGW = energies.EnergyIGWPlus + energies.EnergyIGWMinus;
energies.EnergyIO = energies.EnergyIOBaroclinic + energies.EnergyIOBarotropic;
energies.EnergyGeostrophic = energies.EnergyGeostrophicBaroclinic + energies.EnergyGeostrophicBarotropic;

energies.EnergyTotal = energies.EnergyIGW + energies.EnergyIO + energies.EnergyGeostrophic;

% energies.EnergyDepthIntegrated = ncread(file,'EnergyDepthIntegrated');

netcdfTools = WaveVortexModelNetCDFTools(file);
wvm = netcdfTools.InitializeWaveVortexModelFromNetCDFFile();

t=ncread(file,'t');
% t = t(1:21);
nT = length(t);
dWaveEnergy = zeros(size(t)); dGeostrophicEnergy= zeros(size(t)); dIOEnergy = zeros(size(t));
stride = 1;
for iTime=1:stride:nT
    currentTime = netcdfTools.SetWaveModelToIndex(iTime);
    [Ep,Em,E0] = wvm.EnergyFluxAtTime(currentTime,wvm.Ap,wvm.Am,wvm.A0);
     dIOEnergy(iTime) = 2*sum(Em(1,1,:));
     Ep(1,1,:) = 0; Em(1,1,:) = 0;
     dWaveEnergy(iTime) = sum(Ep(:)) + sum(Em(:));
     dGeostrophicEnergy(iTime) = sum(E0(:));
end

dWaveEnergyActual = diff(energies.EnergyIGW)./diff(energies.t);
dGeostrophicEnergyActual = diff(energies.EnergyGeostrophic)./diff(energies.t);
dEnergyIOActual = diff(energies.EnergyIO)./diff(energies.t);

figure
plot(t(1:stride:nT)/86400,dIOEnergy(1:stride:nT)), hold on
plot(t(1:stride:nT)/86400,dWaveEnergy(1:stride:nT))
plot(t(1:stride:nT)/86400,dGeostrophicEnergy(1:stride:nT))
plot(t(2:end)/86400,dWaveEnergyActual,"Color",'black','LineWidth',2)
plot(t(2:end)/86400,dGeostrophicEnergyActual,"Color",'black','LineWidth',2)
plot(t(2:end)/86400,dEnergyIOActual,"Color",'black','LineWidth',2)
legend('IO','wave','geostrophic', 'wave actual', 'geostrophic actual')
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Internal gravity wave energy
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% total wave energy fluxes
waveEnergyPlus = cumtrapz(t,sum(sum(dWaveEnergyPlusTKJ,3),2));
waveEnergyMinus = cumtrapz(t,sum(sum(dWaveEnergyMinusTKJ,3),2));
dWaveEnergyPredicted = diff(waveEnergyPlus+waveEnergyMinus)./diff(t);

waveEnergyPlus_IO_g = cumtrapz(t,sum(sum(dWaveEnergyPlusTKJ_IO_g,3),2));
waveEnergyMinus_IO_g = cumtrapz(t,sum(sum(dWaveEnergyMinusTKJ_IO_g,3),2));
dWaveEnergyPredicted_IO_g = diff(waveEnergyPlus_IO_g+waveEnergyMinus_IO_g)./diff(t);

% measured wave energy flux
dWaveEnergyActual = diff(energies.EnergyIGW)./diff(energies.t);

FigureSize = [50 50 figure_width_2col+8 220*scaleFactor];
fig1 = figure('Units', 'points', 'Position', FigureSize);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];

sp2 = subplot(2,1,1);

h1=plot(energies.t(2:end)/86400,1e8*dWaveEnergyActual,'k','LineWidth',2); hold on
h2=plot(t(2:end)/86400,1e8*dWaveEnergyPredicted,':k','LineWidth',2);
ax = gca; ax.ColorOrderIndex = 1;
% plot(t(2:end)/86400,dWaveEnergyPredicted_IO_zeta)
h3=plot(t(2:end)/86400,1e8*dWaveEnergyPredicted_IO_g,'LineWidth',2);
% plot(t(2:end)/86400,zeros(size(t(2:end))),'k','LineWidth',1);
xlim([min(t(2:end)/86400) max(t(2:end)/86400)])
xlabel('time (days)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('(m^3/s^2)/s \cdot 10^{-8}', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
legend([h1 h2 h3],'$\frac{d}{dt} E_\textrm{igw} $','total igw flux','$\textrm{io}\nabla\textrm{g}$ flux','Location','south','Interpreter','latex','FontSize',figure_legend_size)
set( gca, 'FontSize', figure_axis_tick_size);
% sp2.YAxisLocation = 'right';
% ylim([0 3e-8])
ylim([0 3])

text(2,2.5,'(a)','FontSize',  figure_axis_label_size, 'FontName', figure_font)

% print('energy-flux-internal-gravity-wave.eps','-depsc2');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Inertial energy
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inertialEnergyPredicted  = cumtrapz(t,sum(dIOEnergyTJ,2));
dInertialEnergyPredicted = diff(inertialEnergyPredicted)./diff(t);

inertialEnergy_wave_g = cumtrapz(t,sum(dIOEnergyTJ_wave_g,2));
dInertialEnergyPredicted_wave_g= diff(inertialEnergy_wave_g)./diff(t);

inertialIOEnergy_wave_wave = cumtrapz(t,sum(dIOEnergyTJ_wave_wave,2));
dInertialIOEnergyPredicted_wave_wave= diff(inertialIOEnergy_wave_wave)./diff(t);

dInertialEnergyActual = diff(energies.EnergyIOBaroclinic)./diff(energies.t);

sp4 = subplot(2,1,2);

h1 = plot(energies.t(2:end)/86400,1e8*dInertialEnergyActual,'k','LineWidth',2); hold on
h2 = plot(t(2:end)/86400,1e8*dInertialEnergyPredicted,':k','LineWidth',2);
ax = gca; ax.ColorOrderIndex = 3;
h3 = plot(t(2:end)/86400,1e8*dInertialEnergyPredicted_wave_g,'LineWidth',2);
h4 = plot(t(2:end)/86400,1e8*dInertialIOEnergyPredicted_wave_wave,'LineWidth',2);
% plot(t(2:end)/86400,zeros(size(t(2:end))),'k','LineWidth',1);
xlim([min(t(2:end)/86400) max(t(2:end)/86400)])
xlabel('time (days)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
% ylabel('energy flux (m^3/s^2)/s', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
legend([h1 h2 h3 h4],'$\frac{d}{dt} E_\textrm{io} $','total io flux','$\textrm{igw}\nabla\textrm{g}$ flux','$\textrm{igw}\nabla\textrm{igw}$ flux','Interpreter','latex','Location','north','FontSize',figure_legend_size)
set( gca, 'FontSize', figure_axis_tick_size);
% sp4.YAxisLocation = 'right';
% ylim([-3e-8 0])
ylim([-3 0])
yticks([-3 -2 -1])
ylabel('energy flux', 'FontSize', figure_axis_label_size, 'FontName', figure_font)

text(2,-2.5,'(b)','FontSize',  figure_axis_label_size, 'FontName', figure_font)



% print('energy-flux-inertial.eps','-depsc2');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % Geostrophic energy
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% geostrophicEnergyPredicted  = cumtrapz(t,sum(sum(dGeostrophicEnergyTKJ,3),2));
% dGeostrophicEnergyPredicted = diff(geostrophicEnergyPredicted)./diff(t);
% 
% dGeostrophicEnergyActual = diff(energies.EnergyGeostrophicBaroclinic)./diff(energies.t);
% 
% sp3 = subplot(2,2,3);
% 
% plot(energies.t(2:end)/86400,dGeostrophicEnergyActual,'k','LineWidth',2), hold on
% plot(t(2:end)/86400,dGeostrophicEnergyPredicted,':k','LineWidth',2)
% xlim([min(t(2:end)/86400) max(t(2:end)/86400)])
% xlabel('time (days)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
% ylabel('energy flux (m^3/s^2)/s', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
% legend('$\frac{d}{dt} E_\textrm{g} $','total flux','Location','southwest','FontSize',figure_legend_size,'Interpreter','latex')
% set( gca, 'FontSize', figure_axis_tick_size);
% ylim([-3e-8 3e-8])
% 
% text(2,2.5e-8,'(c)','FontSize',  figure_axis_label_size, 'FontName', figure_font)
% 
% % print('energy-flux-geostrophic.eps','-depsc2');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % Total energy
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% sp1 = subplot(2,2,1);
% 
% plot(energies.t(2:end)/86400,dWaveEnergyActual+dInertialEnergyActual+dGeostrophicEnergyActual,'k','LineWidth',2), hold on
% plot(t(2:end)/86400,dWaveEnergyPredicted+dInertialEnergyPredicted+dGeostrophicEnergyPredicted,':k','LineWidth',2)
% xlim([min(t(2:end)/86400) max(t(2:end)/86400)])
% xlabel('time (days)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
% ylabel('energy flux (m^3/s^2)/s', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
% legend('$\frac{d}{dt} E_\textrm{total} $','total flux','Location','southwest','FontSize',figure_legend_size,'Interpreter','latex')
% set( gca, 'FontSize', figure_axis_tick_size);
% ylim([-3e-8 3e-8])
% 
% text(2,2.5e-8,'(a)','FontSize',  figure_axis_label_size, 'FontName', figure_font)
% 
% packfig(2,2)
% sp2.YLabel.String = 'energy flux (m^3/s^2)/s';
% sp2.YTickLabelMode = 'auto';
% sp4.YLabel.String = 'energy flux (m^3/s^2)/s';
% sp4.YTickLabelMode = 'auto';

packfig(2,1)
tightfig

sp2.YLabel.Position(1) = -1.0;
sp2.YLabel.Position(2) = 1.0;
sp4.YLabel.Position(1) = -1.0;
sp4.YLabel.Position(2) = -1.0;

print('energy-flux-four-panel.eps','-depsc2');
