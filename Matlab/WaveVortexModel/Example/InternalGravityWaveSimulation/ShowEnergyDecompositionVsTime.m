file = '/Volumes/MoreStorage/Data/InternalWaveSimulation/igw-simulation.nc';

t = ncread(file,'t');

internalWaveEnergyPlus = ncread(file,'EnergyIGWPlus');
internalWaveEnergyMinus = ncread(file,'EnergyIGWMinus');
geostrophicEnergyBaroclinic = ncread(file,'EnergyGeostrophicBaroclinic');
geostrophicEnergyBarotropic = ncread(file,'EnergyGeostrophicBarotropic');
inertialEnergyBaroclinic = ncread(file,'EnergyIOBaroclinic');
inertialEnergyBarotropic = ncread(file,'EnergyIOBarotropic');

WaveEnergyTotal = internalWaveEnergyPlus + internalWaveEnergyMinus;
GeostrophicTotal = geostrophicEnergyBaroclinic + geostrophicEnergyBarotropic;
InertialTotal = inertialEnergyBaroclinic + inertialEnergyBarotropic;

SpectralTotal = WaveEnergyTotal+GeostrophicTotal+InertialTotal;

FigureSize = [50 50 180 230];
fig1 = figure('Units', 'points', 'Position', FigureSize);
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];

plot(t/86400,SpectralTotal/SpectralTotal(1),'LineWidth', 2), hold on
plot(t/86400,InertialTotal/SpectralTotal(1),'LineWidth', 2)
plot(t/86400,WaveEnergyTotal/SpectralTotal(1),'LineWidth', 2)
plot(t/86400,GeostrophicTotal/SpectralTotal(1),'LineWidth', 2)

xlim([0 max(t/86400)])
% xlabel('time (days)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
% ylabel('total energy fraction', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
% set( gca, 'FontSize', figure_axis_tick_size);

legend('total','inertial','wave','geostrophic', 'Location', 'Northeast')

figure
plot(t/86400,GeostrophicTotal/SpectralTotal(1),'LineWidth', 2)
xlabel('time (days)')
ylabel('geostrophic energy')

print('-depsc2', 'geostrophic-energy-vs-time.eps')