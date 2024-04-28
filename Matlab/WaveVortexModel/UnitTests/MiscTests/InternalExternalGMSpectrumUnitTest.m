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

t = 3600;

wvm = WaveVortexModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);
%wavemodel = InternalWaveModelExponentialStratification([Lx, Ly, Lz], [Nx, Ny, Nz], [N0 1300], linspace(-Lz,0,Nz)', latitude);
wvm.initWithGMSpectrum(1.0)
% wavemodel.initWithWaveModes(0,0,1,0,1.0,1);
% [u_gridded,v_gridded] = wavemodel.VelocityFieldAtTime(t);
[u_gridded,v_gridded,w_gridded,rho_prime_gridded] = wvm.VariableFieldsAtTime(t,'u','v','w','rho_prime');

[omega, alpha, k, l, mode, phi, A] = wvm.waveModesFromWaveCoefficients();
L_gm = 1.3e3; % thermocline exponential scale, meters
invT_gm = 5.2e-3; % reference buoyancy frequency, radians/seconds
E_gm = 6.3e-5; % non-dimensional energy parameter
E = L_gm*L_gm*L_gm*invT_gm*invT_gm*E_gm;

fprintf('The wave components sum to %.2f%% GM\n', 100*(wvm.waveEnergy + wvm.inertialEnergy)/E);

% Now reset the gridded field to zero.
wvm.removeAllWaves();

% wavemodel.setExternalWavesWithFrequencies(omega,alpha,J,phi_plus,A_plus,Normalization.kConstant);
wvm.setExternalWavesWithFrequencies(omega, alpha, mode, phi, A,Normalization.kConstant);

[u_ext,v_ext,w_ext,rho_prime_ext] = wvm.VariableFieldsAtTime(t,'u','v','w','rho_prime');

error = @(u,u_unit) max( [max(max(max(abs(u-u_unit)/max( [max(max(max( u ))), 1e-15] )))), 1e-15]);
u_error = error(u_ext,u_gridded);
v_error = error(v_ext,v_gridded);
w_error = error(w_ext,w_gridded);
rho_error = error(rho_prime_ext,rho_prime_gridded);
fprintf('The gridded solution for (u,v) matches the external solution to 1 part in (10^%d, 10^%d) at time t=%d\n', round((log10(u_error))), round((log10(v_error))),t);
fprintf('The gridded solution for (w,rho_prime) matches the external solution to 1 part in (10^%d, 10^%d) at time t=%d\n', round((log10(w_error))), round((log10(rho_error))),t);
