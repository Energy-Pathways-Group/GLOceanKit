%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute Ertel PV
%
% This script uses the InternalWaveModel to create a Garrett-Munk spectrum
% and then compute the Ertel PV
%
% Jeffrey J. Early
% jeffrey@jeffreyearly.com
%
% December 10th, 2018      Version 1.0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aspectRatio = 4;

Lx = 25e3;
Ly = aspectRatio*Lx;
Lz = 1300;

Nx = 4;
Ny = aspectRatio*Nx;
Nz = 9; % 2^n + 1 grid points, to match the Winters model, but 2^n ok too.

latitude = 31;
N0 = 5.2e-3; % Choose your stratification 7.6001e-04

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wavemodel = InternalWaveModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);
% wavemodel.InitializeWithGMSpectrum(1.0);

wavemodel.InitializeWithPlaneWave(0,2,1,0.20,1);

[u,v,w,rho_prime,eta] = wavemodel.VariableFieldsAtTime(0,'u', 'v', 'w', 'rho_prime','zeta');

x = wavemodel.x;
y = wavemodel.y;
z = wavemodel.z;
f0 = wavemodel.f0;

zeta_x = DiffFourier(y,w) - DiffCosine(z,v); % w_y - v_z
zeta_y = DiffCosine(z,u) - DiffFourier(x,w); % u_z - w_x
zeta_z = DiffFourier(x,v) - DiffFourier(y,u); % v_x - u_y

% scaling b so that it is meters (isopycnal height), with a sign difference
b = -(wavemodel.g/wavemodel.rho0)*rho_prime/N0/N0;

% The derivatives are unitless
b_x = DiffFourier(x,b);
b_y = DiffFourier(y,b);
b_z = DiffSine(z,b);

% this should be 1, if we've done this correctly
bbar_z = wavemodel.N2AtDepth(wavemodel.Z)/N0/N0;

PV_x = zeta_x .* b_x;
PV_y = zeta_y .* b_y;
PV_z = (zeta_z + f0) .* b_z + zeta_z .* bbar_z;

PV_linear_1 = zeta_z .* bbar_z;
PV_linear_2 = f0 * b_z;

PV_linear = zeta_z .* bbar_z + f0 * b_z;
PV_z_nonlinear = zeta_z .* b_z;

PV = PV_x + PV_y + PV_z;
