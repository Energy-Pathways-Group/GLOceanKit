file = '/Users/jearly/Data/InternalWaveSimulation/igw-simulation.nc';

netcdfTools = WaveVortexModelNetCDFTools(file);
wvm = netcdfTools.InitializeWaveVortexModelFromNetCDFFile();

t = netcdfTools.SetWaveModelToIndex(2);

[Ep_igw_igw,Em_igw_igw,E0_igw_igw] = wvm.EnergyFluxForFlowConstituentsAtTime(t,wvm.Ap,wvm.Am,wvm.A0,FlowConstituents('internalGravityWave'),FlowConstituents('internalGravityWave'));
[Ep_io_igw,Em_io_igw,E0_io_igw] = wvm.EnergyFluxForFlowConstituentsAtTime(t,wvm.Ap,wvm.Am,wvm.A0,FlowConstituents('inertial'),FlowConstituents('internalGravityWave'));
