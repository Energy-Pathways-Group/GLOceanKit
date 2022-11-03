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

outputfile = '/Users/jearly/Data/InternalWaveSimulation/igw-simulation.nc';
% outputfile = '/Volumes/MoreStorage/Data/InternalWaveSimulation/igw-simulation-constit-fluxes.nc';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wvt = WVTransformConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], N0,latitude=latitude);
wvt.initWithGMSpectrum(1.0)
wvt.summarizeEnergyContent;

outputVar = WVVariableAnnotation('E0_igw_igw',{},'m^3/s^2', 'flux into E0, from wave grad-wave');
f = @(wvt) sum(sum(sum(wvt.energyFluxForFlowConstituents(WVFlowConstituent('wave'),WVFlowConstituent('wave')))));
wvt.addOperation(WVOperation('E0_igw_igw',outputVar,f));

outputVar = WVVariableAnnotation('E0_igw_io',{},'m^3/s^2', 'flux into E0, from wave grad-inertial');
f = @(wvt) sum(sum(sum(wvt.energyFluxForFlowConstituents(WVFlowConstituent('wave'),WVFlowConstituent('inertial')))));
wvt.addOperation(WVOperation('E0_igw_io',outputVar,f));

outputVar = WVVariableAnnotation('E0_igw_g',{},'m^3/s^2', 'flux into E0, from wave grad-geostrophic');
f = @(wvt) sum(sum(sum(wvt.energyFluxForFlowConstituents(WVFlowConstituent('wave'),WVFlowConstituent('geostrophic')))));
wvt.addOperation(WVOperation('E0_igw_g',outputVar,f));

outputVar = WVVariableAnnotation('E0_io_igw',{},'m^3/s^2', 'flux into E0, from inertial grad-wave');
f = @(wvt) sum(sum(sum(wvt.energyFluxForFlowConstituents(WVFlowConstituent('inertial'),WVFlowConstituent('wave')))));
wvt.addOperation(WVOperation('E0_io_igw',outputVar,f));

outputVar = WVVariableAnnotation('E0_io_g',{},'m^3/s^2', 'flux into E0, from inertial grad-geostrophic');
f = @(wvt) sum(sum(sum(wvt.energyFluxForFlowConstituents(WVFlowConstituent('inertial'),WVFlowConstituent('geostrophic')))));
wvt.addOperation(WVOperation('E0_io_g',outputVar,f));

outputVar = WVVariableAnnotation('E0_g_io',{},'m^3/s^2', 'flux into E0, from geostrophic grad-wave');
f = @(wvt) sum(sum(sum(wvt.energyFluxForFlowConstituents(WVFlowConstituent('geostrophic'),WVFlowConstituent('wave')))));
wvt.addOperation(WVOperation('E0_g_io',outputVar,f));

model = WVModel(wvt,nonlinearFlux=BoussinesqConstantN(wvt,shouldAntialias=1));
model.setupIntegrator(timeStepConstraint="advective", outputInterval=wvt.inertialPeriod/10);
model.createNetCDFFileForModelOutput(outputfile,shouldOverwriteExisting=1);
model.addNetCDFOutputVariables('E0_igw_igw','E0_igw_io','E0_igw_g','E0_io_igw','E0_io_g','E0_g_io');
model.addNetCDFOutputVariables('geostrophicEnergyBaroclinic','geostrophicEnergyBarotropic')
model.integrateToTime(3*wvt.inertialPeriod);

ncfile = model.ncfile;
t = ncfile.readVariables('t');

figure
plot(t/inertialPeriod,E0_all_all,'LineWidth', 4,'Color','k'), hold on
plot(t/inertialPeriod,E0_igw_igw,'LineWidth', 2)
plot(t/inertialPeriod,E0_igw_io,'LineWidth', 2)
plot(t/inertialPeriod,E0_igw_g,'LineWidth', 2)
plot(t/inertialPeriod,E0_io_igw,'LineWidth', 2)
plot(t/inertialPeriod,E0_io_g,'LineWidth', 2)
plot(t/inertialPeriod,E0_g_igw,'LineWidth', 2)
legend('total','igw-igw','igw-io','igw-g','io-igw','io-g','g-igw')
xlabel('time (inertial periods)')
ylabel('energy flux into the geostrophic field')
% print('EnergyFluxIntoPV.eps','-depsc2')

[geostrophicEnergyBaroclinic,geostrophicEnergyBarotropic] = ncfile.readVariables('geostrophicEnergyBaroclinic','geostrophicEnergyBarotropic');
GeostrophicTotal = geostrophicEnergyBaroclinic + geostrophicEnergyBarotropic;
figure
plot(t(2:end)/inertialPeriod,diff(GeostrophicTotal)./diff(t))
xlabel('time (inertial periods)')
ylabel('change in geostrophic energy')
% print('ChangeInPV.eps','-depsc2')


