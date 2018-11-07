%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GMExponentialUnitTest
%
% This script uses the InternalWaveModel to create, and validate, a
% Garrett-Munk spectrum in a linear internal wave field with exponential
% stratification.
%
% Jeffrey J. Early
% jeffrey@jeffreyearly.com
%
% November 5th, 2018      Version 1.0



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aspectRatio = 4;

Lx = 100e3;
Ly = aspectRatio*100e3;
Lz = 5000;

N = 32;
Nx = N;
Ny = aspectRatio*N;
Nz = N+1;

latitude = 31;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wavemodel = InternalWaveModelExponentialStratification([Lx, Ly, Lz], [Nx, Ny, Nz], [5.2e-3 1300], linspace(-Lz,0,Nz), latitude);
% wavemodel.FillOutWaveSpectrum();
wavemodel.InitializeWithGMSpectrum(1.0);

[u,v,w] = wavemodel.VelocityFieldAtTime(0);
