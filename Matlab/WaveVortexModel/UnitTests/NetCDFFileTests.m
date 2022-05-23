%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 64;

Lx = 50e3;
Ly = Lx;
Lz = 1300;

Nx = N;
Ny = N;
Nz = N+1; % 2^n + 1 grid points, to match the Winters model, but 2^n ok too.

latitude = 25;
N0 = 5.2e-3; % Choose your stratification 7.6001e-04

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0
    wvt = WaveVortexTransformConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);
else
    rho0 = 1025; g = 9.81;
    rho = @(z) -(N0*N0*rho0/g)*z + rho0;
    N2Function = @(z) N0*N0*ones(size(z));
    dLnN2Function = @(z) zeros(size(z));
    wvt = WaveVortexTransformHydrostatic([Lx, Ly, Lz], [Nx, Ny, Nz-1], latitude, rho,'N2func', N2Function, 'dLnN2func',dLnN2Function);
end

U = .2;
period = wvt.InitializeWithPlaneWave(10,0,1,U,1);  

wvt.WriteToFile('test.nc','Ap','Am','A0');

wvt2 = WaveVortexTransform.InitFromNetCDFFile('test.nc',243);