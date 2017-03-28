%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% InternalExternalGMSpectrumUnitTest
%
% This script tests to see if the gridded IW solutions are identical to the
% external IW solutions with the same amplitude and phases.
%
% Jeffrey J. Early
% jeffrey@jeffreyearly.com
%
% March 27th, 2017      Version 1.0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aspectRatio = 1;

Lx = 100e3;
Ly = aspectRatio*100e3;
Lz = 5000;

Nx = 4;
Ny = aspectRatio*4;
Nz = 5; % 2^n + 1 grid points, to match the Winters model, but 2^n ok too.

latitude = 31;
N0 = 5.2e-3; % Choose your stratification 7.6001e-04

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = 0;

wavemodel = InternalWaveModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);
% wavemodel.InitializeWithGMSpectrum(1.0)
wavemodel.InitializeWithPlaneWave(0,0,1,1.0,1);
[u_gridded,v_gridded] = wavemodel.VelocityFieldAtTime(t);

[omega, alpha, mode, phi, A] = wavemodel.WaveCoefficientsFromGriddedWaves();
L_gm = 1.3e3; % thermocline exponential scale, meters
invT_gm = 5.2e-3; % reference buoyancy frequency, radians/seconds
E_gm = 6.3e-5; % non-dimensional energy parameter
E = L_gm*L_gm*L_gm*invT_gm*invT_gm*E_gm;
E = E*(wavemodel.Lz/L_gm);

fprintf('The wave components sum to %.2f%% GM\n', 100*sum(A.*A)/E);

% Now reset the gridded field to zero.
wavemodel.InitializeWithPlaneWave(1,1,1,0.0,1);

% wavemodel.SetExternalWavesWithFrequencies(omega,alpha,J,phi_plus,A_plus,'energyDensity');
wavemodel.SetExternalWavesWithFrequencies(omega, alpha, mode, phi, A,'energyDensity');

[u_ext,v_ext] = wavemodel.VelocityFieldAtTime(t);

error = @(u,u_unit) max( [max(max(max(abs(u-u_unit)/max( [max(max(max( u ))), 1e-15] )))), 1e-15]);
u_error = error(u_ext,u_gridded);
v_error = error(v_ext,v_gridded);
fprintf('The gridded solution for (u,v) matches the external solution to 1 part in (10^%d, 10^%d) at time t=%d\n', round((log10(u_error))), round((log10(v_error))),t);
