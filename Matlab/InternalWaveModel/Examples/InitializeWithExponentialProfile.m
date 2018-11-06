%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Lx = 100e3;
Ly = 100e3;
Lz = 5000;

Nx = 64;
Ny = 64;
Nz = 65; % 2^n + 1 grid points, to match the Winters model, but 2^n ok too.

latitude = 31;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[rhoFunction, N2Function, zIn] = InternalModes.StratificationProfileWithName('exponential');

z = linspace(-Lz,0,Nz)';
% wavemodel = InternalWaveModelArbitraryStratification([Lx, Ly, Lz], [Nx, Ny, Nz], rhoFunction, z, [], latitude);
wavemodel = InternalWaveModelExponentialStratification([Lx, Ly, Lz], [Nx, Ny, Nz], [5.2e-3 1300], z, [], latitude);
wavemodel.InitializeWithGMSpectrum(1.0);