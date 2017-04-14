Lx = 150e3;
Ly = 150e3;
Lz = 5000;

Nx = 64;
Ny = 64;
Nz = 64;

latitude = 31;
N0 = 5.2e-3; % Choose your stratification 7.6001e-04

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wavemodel = InternalWaveModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);
wavemodel.FillOutWaveSpectrum();
wavemodel.InitializeWithGMSpectrum(1.0);

timeStep = 1; % in seconds
maxTime = 100;
t = (0:timeStep:maxTime)';

tic
for iTime=1:length(t)
    [u,v,w]=wavemodel.VelocityFieldAtTime(t(iTime));
    [w,zeta] = wavemodel.VerticalFieldsAtTime(t(iTime));
end
toc