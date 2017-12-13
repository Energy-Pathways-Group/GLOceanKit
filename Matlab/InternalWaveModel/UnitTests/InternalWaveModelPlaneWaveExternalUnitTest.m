%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% InternalWaveModelPlaneWaveUnitTest
%
% This script uses the InternalWaveModel to create, and validate, a single
% internal wave.
%
% Jeffrey J. Early
% jeffrey@jeffreyearly.com
%
% March 25th, 2016      Version 1.0
% March 30th, 2016      Version 1.1
% November 17th, 2016   Version 1.2


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lx = 15e3;
Ly = 15e3;
Lz = 5000;

Nx = 32;
Ny = 32;
Nz = 32;

latitude = 31;
N0 = 5.2e-3; % Choose your stratification 7.6001e-04

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wavemodel = InternalWaveModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);

% rho0 = 1025; g = 9.81;
% rho = @(z) -(N0*N0*rho0/g)*z + rho0;
% z = (Lz/Nz)*(0:Nz-1)' - Lz;
% wavemodel = InternalWaveModelArbitraryStratification([Lx, Ly, Lz], [Nx, Ny, Nz], rho, z, Nz, latitude);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create a single plane-wave with the model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j0 = 20; % j=1..nModes, where 1 indicates the 1st baroclinic mode
U = 0.01; % m/s
sign = 1;
phi = 0;



omega = 2*wavemodel.f0;
alpha = 0;
k = wavemodel.SetExternalWavesWithFrequencies(omega,alpha,j0,phi,U,Normalization.uMax);
k0 = k*cos(alpha);
l0 = k*sin(alpha);

t = 4*86400;
[u,v,w] = wavemodel.VelocityFieldAtTime(t);
zeta = wavemodel.IsopycnalDisplacementFieldAtTime(t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create a single plane-wave with the known analytical solution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f0 = wavemodel.f0;
k = wavemodel.k;
l = wavemodel.l;
h = wavemodel.h;
x = wavemodel.x;
y = wavemodel.y;
X = wavemodel.X;
Y = wavemodel.Y;
Z = wavemodel.Z;

m = j0*pi/Lz;

alpha=atan2(l0,k0);
K = sqrt( k0^2 + l0^2);
theta = k0*X + l0*Y + omega*t;
u_unit = U*(cos(alpha)*cos( theta ) + (f0/omega)*sin(alpha)*sin( theta )).*cos(m*Z);
v_unit = U*(sin(alpha)*cos( theta ) - (f0/omega)*cos(alpha)*sin( theta )).*cos(m*Z);
w_unit = (U*K/m) * sin(theta) .* sin(m*Z);
zeta_unit = -(U*K/m/omega) * cos(theta) .* sin(m*Z);

% This is a measure of relative error, bounded by 1e-15
error = @(u,u_unit) max( [max(max(max(abs(u-u_unit)/max( [max(max(max( u ))), 1e-15] )))), 1e-15]);

% Compute the relative error
u_error = error(u,u_unit);
v_error = error(v,v_unit);
w_error = error(w,w_unit);
zeta_error = error(zeta,zeta_unit);

max_error = max([round((log10(u_error)))  round((log10(v_error))) round((log10(w_error))) round((log10(zeta_error)))]);

fprintf('Testing IW mode (k0,l0,j0)=(%d,%d,%d):\n',k0,l0,j0);
fprintf('The model solution for (u,v) matches the analytical solution to 1 part in (10^%d, 10^%d) at time t=%d\n', round((log10(u_error))), round((log10(v_error))),t);
fprintf('The model solution for (w,zeta) matches the analytical solution to 1 part in (10^%d, 10^%d) at time t=%d\n', round((log10(w_error))), round((log10(zeta_error))),t);
