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

outputfile = '/Users/jearly/Data/InternalWaveSimulation/igw-simulation.nc';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wvt = WVTransformConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], N0,latitude=latitude);
wvt.initWithGMSpectrum(1.0)
wvt.removeEnergyFromAliasedModes();
wvt.summarizeEnergyContent;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Setup a netcdf file for the output
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model = WVModel(wvt,nonlinearFlux=BoussinesqConstantN(wvt,shouldAntialias=1));
model.setupIntegrator(timeStepConstraint="advective", outputInterval=wvt.inertialPeriod/10);
% model.createNetCDFFileForModelOutput(outputfile,shouldOverwriteExisting=1);
model.integrateToTime(3*wvt.inertialPeriod);

% ncfile = model.ncfile;
% t = ncfile.readVariables('t');

% deltaT = 175;
% totalOuterLoop = round(2*wvt.inertialPeriod/deltaT/5);
% 
% integrator2 = ArrayIntegrator(@(t,y0) nonlinearFluxAtTime(wvt,t,y0),{wvt.Ap,wvt.Am,wvt.A0},deltaT);
% wvt.summarizeEnergyContent;
% for j=1:totalOuterLoop
% for i=1:5
%     integrator2.IncrementForward();
%     wvt.Ap = integrator2.currentY{1};
%     wvt.Am = integrator2.currentY{2};
%     wvt.A0 = integrator2.currentY{3};
% end
% wvm.summarizeEnergyContent;
% wvt.summarizeEnergyContent;
% end
% 
% 
% function Farray = nonlinearFluxAtTime(wvt,t,y0)
% wvt.t = t;
% [Fp,Fm,F0] = wvt.nonlinearFlux;
% Farray = {Fp,Fm,F0};
% end