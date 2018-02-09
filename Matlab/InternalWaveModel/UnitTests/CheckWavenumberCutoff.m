%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CheckWavenumberCutoff
%
% Tests the API for limiting the the horizontal wavenumber used to
% initialize waves.
%
% Jeffrey J. Early
% jeffrey@jeffreyearly.com
%
% February 8, 2018      Version 1.0


N = 16;
aspectRatio = 1;

Lx = 100e3;
Ly = aspectRatio*100e3;
Lz = 5000;

Nx = N;
Ny = aspectRatio*N;
Nz = N+1; % 2^n + 1 grid points, to match the Winters model, but 2^n ok too.

latitude = 31;
N0 = 5.2e-3; % Choose your stratification 7.6001e-04

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wavemodel = InternalWaveModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);
wavemodel.InitializeWithGMSpectrum(1.0)

[omega, alpha, k, l, mode, phi, A] = wavemodel.WaveCoefficientsFromGriddedWaves();
K = sqrt( k.*k + l.*l );
fprintf('min wavenumber: %f, max wavenumber: %f\n', min(K), max(K));

delta = Lx/N;
maxK = 2/delta;
wavemodel.InitializeWithGMSpectrum(1.0,'maxK',maxK);
[omega, alpha, k, l, mode, phi, A] = wavemodel.WaveCoefficientsFromGriddedWaves();
K = sqrt( k.*k + l.*l );
fprintf('min wavenumber: %f, max wavenumber: %f\n', min(K), max(K));