%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CustomSpectrumTest
%
% This script tests the API to set a custom wave spectrum
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

wavemodel = InternalWaveModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create a customized spectral function, let's just use GM
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GM Parameters
j_star = 3;
L_gm = 1.3e3; % thermocline exponential scale, meters
invT_gm = 5.2e-3; % reference buoyancy frequency, radians/seconds
E_gm = 6.3e-5; % non-dimensional energy parameter
E = L_gm*L_gm*L_gm*invT_gm*invT_gm*E_gm;

% Compute the proper vertical function normalization
H = (j_star+(1:1024)).^(-5/2);
H_norm = 1/sum(H);

% Do the same for the frequency function.
f0 = wavemodel.f0;
B_norm = 1/atan(sqrt(wavemodel.Nmax*wavemodel.Nmax/(f0*f0)-1));

% This function tells you how much energy you need between two
% frequencies for a given vertical mode.
GM2D_int = @(omega0,omega1,j) E*H_norm*B_norm*((j+j_star).^(-5/2))*(atan(f0/sqrt(omega0*omega0-f0*f0)) - atan(f0/sqrt(omega1*omega1-f0*f0)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Now initialize with that spectrum
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wavemodel.InitializeWithSpectralFunction(GM2D_int);

% Now everything should be initialized
[u_gridded,v_gridded,w_gridded,rho_prime_gridded] = wavemodel.VariableFieldsAtTime(t,'u','v','w','rho_prime');

% And we can be cute and pull out the wave amplitude and sum them.
[omega, alpha, mode, phi, A] = wavemodel.WaveCoefficientsFromGriddedWaves();
fprintf('The wave components sum to %.2f%% GM\n', 100*sum(A.*A/2)/E);