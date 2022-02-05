%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GMSpectrum
%
% Time steps a internal wave spectrum using the WaveVortexModel.
%
% Jeffrey J. Early
% jeffrey@jeffreyearly.com
%
% February 3rd, 2022       Version 1.0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 64;
aspectRatio = 1;

Lx = 100e3;
Ly = aspectRatio*100e3;
Lz = 1300;

Nx = N;
Ny = aspectRatio*N;
Nz = N+1; % 2^n + 1 grid points, to match the Winters model, but 2^n ok too.

latitude = 31;
N0 = 5.2e-3; % Choose your stratification 7.6001e-04
rho0 = 1025;

% Simulation length
inertialPeriod = (2*pi/(2 * 7.2921E-5 * sin( latitude*pi/180 )));
maxTime = 2*inertialPeriod;
outputInterval = inertialPeriod/10;
cfl = .25;

outputfile = '/Users/jearly/Data/InternalWaveSimulation/igw-simulation.nc';
outputfile = '/Volumes/MoreStorage/Data/InternalWaveSimulation/igw-simulation-constit-fluxes.nc';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = 3600;

wvm = WaveVortexModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0,rho0,'shouldAntiAlias',1);
% wvm.shouldAntiAlias = 1;
wvm.InitializeWithGMSpectrum(1.0)
wvm.summarizeEnergyContent;

deltaT = wvm.TimeStepForCFL(cfl,outputInterval);
% stepsPerOutput = round(outputInterval/deltaT);
stepsPerOutput = 1;
maxTime = 100*deltaT;
t = (0:deltaT:maxTime)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Setup a netcdf file for the output
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

netcdfTool = WaveVortexModelNetCDFTools(outputfile);
netcdfTool.CreateNetCDFFileFromModel(wvm,length(t),'double');
netcdfTool.CreateAmplitudeCoefficientVariables();
netcdfTool.CreateEnergeticsVariables();
netcdfTool.CreateEnergeticsKJVariables();

% Save the initial conditions
iTime = 1;
netcdfTool.WriteTimeAtIndex(iTime,t(iTime));
netcdfTool.WriteAmplitudeCoefficientsAtIndex(iTime);
netcdfTool.WriteEnergeticsAtIndex(iTime);
netcdfTool.WriteEnergeticsKJAtIndex(iTime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Time step forward!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

integrator = ArrayIntegrator(@(t,y0) wvm.NonlinearFluxAtTimeArray(t,y0),{wvm.Ap,wvm.Am,wvm.A0},deltaT);

E0_all_all_t= zeros(length(t),1);
E0_igw_igw_t = zeros(length(t),1);
E0_igw_io_t = zeros(length(t),1);
E0_igw_g_t = zeros(length(t),1);
E0_io_igw_t = zeros(length(t),1);
E0_io_g_t = zeros(length(t),1);
E0_g_igw_t = zeros(length(t),1);

iTime = 1;
[Ep,Em,E0] = wvm.EnergyFluxAtTime(t(iTime),wvm.Ap,wvm.Am,wvm.A0);
E0_all_all_t(iTime) = sum(E0(:));
[Ep,Em,E0] = wvm.EnergyFluxForFlowConstituentsAtTime(t(iTime),wvm.Ap,wvm.Am,wvm.A0,FlowConstituents('wave'),FlowConstituents('wave'));
E0_igw_igw_t(iTime) = sum(E0(:));
[Ep,Em,E0] = wvm.EnergyFluxForFlowConstituentsAtTime(t(iTime),wvm.Ap,wvm.Am,wvm.A0,FlowConstituents('wave'),FlowConstituents('inertial'));
E0_igw_io_t(iTime) = sum(E0(:));
[Ep,Em,E0] = wvm.EnergyFluxForFlowConstituentsAtTime(t(iTime),wvm.Ap,wvm.Am,wvm.A0,FlowConstituents('wave'),FlowConstituents('geostrophic'));
E0_igw_g_t(iTime) = sum(E0(:));
[Ep,Em,E0] = wvm.EnergyFluxForFlowConstituentsAtTime(t(iTime),wvm.Ap,wvm.Am,wvm.A0,FlowConstituents('inertial'),FlowConstituents('wave'));
E0_io_igw_t(iTime) = sum(E0(:));
[Ep,Em,E0] = wvm.EnergyFluxForFlowConstituentsAtTime(t(iTime),wvm.Ap,wvm.Am,wvm.A0,FlowConstituents('wave'),FlowConstituents('geostrophic'));
E0_io_g_t(iTime) = sum(E0(:));
[Ep,Em,E0] = wvm.EnergyFluxForFlowConstituentsAtTime(t(iTime),wvm.Ap,wvm.Am,wvm.A0,FlowConstituents('geostrophic'),FlowConstituents('wave'));
E0_g_igw_t(iTime) = sum(E0(:));

% profile on
for iTime=2:length(t)
    if iTime == 2
        startTime = datetime('now');
        fprintf('Starting numerical simulation on %s\n', datestr(startTime));
    end
    if iTime == 3 || mod(iTime,5) == 0
        timePerStep = (datetime('now')-startTime)/(iTime-2);
        timeRemaining = (length(t)-iTime+1)*timePerStep;
        fprintf('\twriting values time step %d of %d to file. Estimated finish time %s (%s from now)\n', iTime, length(t), datestr(datetime('now')+timeRemaining), datestr(timeRemaining, 'HH:MM:SS')) ;
    end

    for i=1:stepsPerOutput
        integrator.IncrementForward();
        wvm.Ap = integrator.currentY{1};
        wvm.Am = integrator.currentY{2};
        wvm.A0 = integrator.currentY{3};
    end

    [Ep,Em,E0] = wvm.EnergyFluxAtTime(t(iTime),wvm.Ap,wvm.Am,wvm.A0);
    E0_all_all_t(iTime) = sum(E0(:));
    [Ep,Em,E0] = wvm.EnergyFluxForFlowConstituentsAtTime(t(iTime),wvm.Ap,wvm.Am,wvm.A0,FlowConstituents('wave'),FlowConstituents('wave'));
    E0_igw_igw_t(iTime) = sum(E0(:));
    [Ep,Em,E0] = wvm.EnergyFluxForFlowConstituentsAtTime(t(iTime),wvm.Ap,wvm.Am,wvm.A0,FlowConstituents('wave'),FlowConstituents('inertial'));
    E0_igw_io_t(iTime) = sum(E0(:));
    [Ep,Em,E0] = wvm.EnergyFluxForFlowConstituentsAtTime(t(iTime),wvm.Ap,wvm.Am,wvm.A0,FlowConstituents('wave'),FlowConstituents('geostrophic'));
    E0_igw_g_t(iTime) = sum(E0(:));
    [Ep,Em,E0] = wvm.EnergyFluxForFlowConstituentsAtTime(t(iTime),wvm.Ap,wvm.Am,wvm.A0,FlowConstituents('inertial'),FlowConstituents('wave'));
    E0_io_igw_t(iTime) = sum(E0(:));
    [Ep,Em,E0] = wvm.EnergyFluxForFlowConstituentsAtTime(t(iTime),wvm.Ap,wvm.Am,wvm.A0,FlowConstituents('wave'),FlowConstituents('geostrophic'));
    E0_io_g_t(iTime) = sum(E0(:));
    [Ep,Em,E0] = wvm.EnergyFluxForFlowConstituentsAtTime(t(iTime),wvm.Ap,wvm.Am,wvm.A0,FlowConstituents('geostrophic'),FlowConstituents('wave'));
    E0_g_igw_t(iTime) = sum(E0(:));

    if mod(iTime,10)==0
        wvm.summarizeEnergyContent();
    end
    netcdfTool.WriteTimeAtIndex(iTime,t(iTime));
    netcdfTool.WriteAmplitudeCoefficientsAtIndex(iTime);
    netcdfTool.WriteEnergeticsAtIndex(iTime);
    netcdfTool.WriteEnergeticsKJAtIndex(iTime);
end

netcdfTool.close();

figure
plot(t/inertialPeriod,E0_all_all_t,'LineWidth', 4,'Color','k'), hold on
plot(t/inertialPeriod,E0_igw_igw_t,'LineWidth', 2)
plot(t/inertialPeriod,E0_igw_io_t,'LineWidth', 2)
plot(t/inertialPeriod,E0_igw_g_t,'LineWidth', 2)
plot(t/inertialPeriod,E0_io_igw_t,'LineWidth', 2)
plot(t/inertialPeriod,E0_io_g_t,'LineWidth', 2)
plot(t/inertialPeriod,E0_g_igw_t,'LineWidth', 2)
legend('total','igw-igw','igw-io','igw-g','io-igw','io-g','g-igw')
xlabel('time (inertial periods)')
ylabel('energy flux into the geostrophic field')
print('EnergyFluxIntoPV.eps','-depsc2')

baroclinicGeostrophicEnergy = ncread(outputfile,'EnergyGeostrophicBaroclinic');
barotropicGeostrophicEnergy = ncread(outputfile,'EnergyGeostrophicBarotropic');
GeostrophicTotal = baroclinicGeostrophicEnergy + barotropicGeostrophicEnergy;
figure
plot(t(2:end)/inertialPeriod,diff(GeostrophicTotal)./diff(t))
xlabel('time (inertial periods)')
ylabel('change in geostrophic energy')
print('ChangeInPV.eps','-depsc2')


