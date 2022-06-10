file = '/Users/jearly/Data/InternalWaveSimulation/igw-simulation.nc';

netcdfTools = WaveVortexModelNetCDFTools(file);
wvm = netcdfTools.InitializeWaveVortexModelFromNetCDFFile();
Nt = netcdfTools.Nt;

E0_igw_igw_t = zeros(Nt,1);
E0_io_igw_t = zeros(Nt,1);

for iTime = 1:netcdfTools.Nt
    t = netcdfTools.SetWaveModelToIndex(iTime);

    [Ep_igw_igw,Em_igw_igw,E0_igw_igw] = wvm.EnergyFluxForFlowConstituentsAtTime(t,wvm.Ap,wvm.Am,wvm.A0,FlowConstituents('internalGravityWave'),FlowConstituents('internalGravityWave'));
    [Ep_io_igw,Em_io_igw,E0_io_igw] = wvm.EnergyFluxForFlowConstituentsAtTime(t,wvm.Ap,wvm.Am,wvm.A0,FlowConstituents('inertial'),FlowConstituents('internalGravityWave'));

    E0_igw_igw_t(iTime) = sum(E0_igw_igw(:));
    E0_io_igw_t(iTime) = sum(E0_io_igw(:));
end

t = ncread(file,'t');
geostrophicEnergyBaroclinic = ncread(file,'EnergyGeostrophicBaroclinic');
geostrophicEnergyBarotropic = ncread(file,'EnergyGeostrophicBarotropic');
GeostrophicTotal = geostrophicEnergyBaroclinic + geostrophicEnergyBarotropic;

figure
plot(t,E0_igw_igw_t,'LineWidth', 2), hold on
plot(t,E0_io_igw_t,'LineWidth', 2), hold on

plot(t(2:end),diff(GeostrophicTotal)./diff(t))