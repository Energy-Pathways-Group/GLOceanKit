%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nz = 64;
Lz = 1300;

latitude = 25;
N0 = 5.2e-3; % Choose your stratification 7.6001e-04

N2 = @(z) N0*N0*ones(size(z));

flow = WVStratifiedFlowHydrostatic(Lz, Nz, N2=N2,latitude=latitude);
wvt.initWithRandomFlow(uvMax=0.05);

wvt.nonlinearFluxOperation.addForcing(WVSpectralVanishingViscosity(wvt,uv_damp=wvt.uvMax))

model = WVModel(wvt);

% wvt.writeToFile('test.nc',shouldOverwriteExisting=1);
% 
% wvt2 = WVTransform.waveVortexTransformFromFile("test.nc");

% wvt.nonlinearFluxOperation = WVHydrostaticFlux(wvt,WVSpectralVanishingViscosity(wvt,uv_damp=0.2));
