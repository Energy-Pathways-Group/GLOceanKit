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

aspectRatio = 1;

Lx = 25e3;
Ly = aspectRatio*Lx;
Lz = 1300;

Nx = 64;
Ny = aspectRatio*Nx;
Nz = 65; % 2^n + 1 grid points, to match the Winters model, but 2^n ok too.

latitude = 31;
N0 = 5.2e-3; % Choose your stratification 7.6001e-04

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wavemodel = InternalWaveModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);
wavemodel.InitializeWithGMSpectrum(1.0);

% U = 0.20;
% k0 = 2;
% l0 = 0;
% j0 = 1;
% period = wavemodel.InitializeWithPlaneWave(k0,l0,j0,U,1);
% omega = 2*pi/period;
% kk = wavemodel.k(k0+1);
% ll = wavemodel.l(l0+1);
% mm = j0*pi/wavemodel.Lz;

[u,v,w,rho_prime,eta] = wavemodel.VariableFieldsAtTime(0,'u', 'v', 'w', 'rho_prime','zeta');

x = wavemodel.x;
y = wavemodel.y;
z = wavemodel.z;
f0 = wavemodel.f0;

zeta_x = DiffFourier(y,w,1,2) - DiffCosine(z,v); % w_y - v_z
zeta_y = DiffCosine(z,u) - DiffFourier(x,w,1,1); % u_z - w_x
zeta_z = DiffFourier(x,v,1,1) - DiffFourier(y,u,1,2); % v_x - u_y

% scaling b so that it is meters (isopycnal height), with a sign difference
b = -(wavemodel.g/wavemodel.rho0)*rho_prime/N0/N0;

% The derivatives are unitless
b_x = DiffFourier(x,b,1,1);
b_y = DiffFourier(y,b,1,2);
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

X = wavemodel.X;
Y = wavemodel.Y;
Z = wavemodel.Z;
cp = omega/kk;
PV_predicted = - f0*(U/cp)^2*(-(sin(kk*X)).^2 + (cos(mm*(Z+Lz))).^2);
