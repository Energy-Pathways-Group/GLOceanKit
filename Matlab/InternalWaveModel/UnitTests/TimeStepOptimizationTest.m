Lx = 150e3;
Ly = 150e3;
Lz = 5000;

Nx = 64;
Ny = 64;
Nz = 64+1;

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

dx = wavemodel.x(2)-wavemodel.x(1);
dy = wavemodel.y(2)-wavemodel.y(1);
N = 10;
x_float = (1:N)*dx;
y_float = (1:N)*dy;
z_float = (0:2)*(-Lz/4);

[x_float,y_float,z_float] = ndgrid(x_float,y_float,z_float);
x_float = reshape(x_float,[],1);
y_float = reshape(y_float,[],1);
z_float = reshape(z_float,[],1);

timeStep = 1; % in seconds
maxTime = 10;
t = (0:timeStep:maxTime)';

tic
for iTime=1:length(t)
%     [u,v,w]=wavemodel.VelocityFieldAtTime(t(iTime));
%     [w,zeta] = wavemodel.VerticalFieldsAtTime(t(iTime));
%     isopycnalDeviation = wavemodel.ZetaAtTimePosition(0,x_float,y_float,z_float);
    [u,v,w] = wavemodel.VelocityAtTimePosition(t(iTime),x_float,y_float,z_float,'exact');
end
toc