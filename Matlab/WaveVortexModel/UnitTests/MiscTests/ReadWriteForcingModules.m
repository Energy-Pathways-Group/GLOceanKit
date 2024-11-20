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
phi=pi/3;

N2 = @(z) N0*N0*ones(size(z));

wvt = WVTransformHydrostatic([Lx, Ly, Lz], [Nx, Ny, Nz], N2=N2,latitude=latitude);
wvt.initWithRandomFlow(uvMax=0.05);

wvt.nonlinearFluxOperation.addForcing(WVSpectralVanishingViscosity(wvt,uv_damp=wvt.uvMax))

model = WVModel(wvt);

% wvt.writeToFile('test.nc',shouldOverwriteExisting=1);
% 
% wvt2 = WVTransform.waveVortexTransformFromFile("test.nc");

% wvt.nonlinearFluxOperation = WVHydrostaticFlux(wvt,WVSpectralVanishingViscosity(wvt,uv_damp=0.2));
