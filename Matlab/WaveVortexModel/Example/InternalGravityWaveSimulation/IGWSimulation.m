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

outputfile = '/Volumes/MoreStorage/Data/InternalWaveSimulation/igw-simulation.nc';

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
stepsPerOutput = round(outputInterval/deltaT);
t = (0:outputInterval:maxTime)';

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

    if mod(iTime,10)==0
        wvm.summarizeEnergyContent();
    end
    netcdfTool.WriteTimeAtIndex(iTime,t(iTime));
    netcdfTool.WriteAmplitudeCoefficientsAtIndex(iTime);
    netcdfTool.WriteEnergeticsAtIndex(iTime);
    netcdfTool.WriteEnergeticsKJAtIndex(iTime);
end

netcdfTool.close();