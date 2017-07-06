Lx = 800e3;
Ly = 800e3;
Lz = 5000;

N = 64;
Nx = N;
Ny = N;
Nz = N+1;

latitude = 31;
N0 = 5.2e-3; % Choose your stratification 7.6001e-04

interpolationMethod = 'spline';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wavemodel = InternalWaveModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);
wavemodel.FillOutWaveSpectrum();
wavemodel.InitializeWithGMSpectrum(Lz/1300);
wavemodel.ShowDiagnostics();

[u,v] = wavemodel.VelocityFieldAtTime(0.0);
U = max(max(max( sqrt(u.*u + v.*v) )));
[w,zeta] = wavemodel.VerticalFieldsAtTime(0.0);
Zeta = max(max(max( zeta )));
fprintf('Max fluid velocity: %.2f cm/s, max isopycnal deviation: %.2f m\n',U*100, Zeta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create floats/drifters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dx = wavemodel.x(2)-wavemodel.x(1);
dy = wavemodel.y(2)-wavemodel.y(1);
nLevels = 5;
N = floor(N/3);
x_float = (0:N-1)*dx;
y_float = (0:N-1)*dy;
z_float = (0:nLevels-1)*(-Lz/(2*(nLevels-1)));

% nudge towards the center of the domain. This isn't necessary, but does
% prevent the spline interpolation from having to worry about the
% boundaries.
x_float = x_float + (max(wavemodel.x) - max(x_float))/2;
y_float = y_float + (max(wavemodel.y) - max(y_float))/2;

[x_float,y_float,z_float] = ndgrid(x_float,y_float,z_float);
x_float = reshape(x_float,[],1);
y_float = reshape(y_float,[],1);
z_float = reshape(z_float,[],1);
nFloats = numel(x_float);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% First try built in, optimized method
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
zIsopycnal1 = wavemodel.PlaceParticlesOnIsopycnal(x_float,y_float,z_float,interpolationMethod,1e-6);
totalTime = toc;
fprintf('Total time %.2f\n',totalTime);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Try generic implementation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
X = wavemodel.X;
Y = wavemodel.Y;
Z = wavemodel.Z;
Rho = wavemodel.DensityAtTime(0);
drhobardz = -wavemodel.N2*wavemodel.rho0/9.81;
zIsopycnal2 = PlaceParticlesOnIsopycnal(x_float,y_float,z_float,X,Y,Z,Rho,drhobardz,interpolationMethod,1e-6);
totalTime = toc;
fprintf('Total time %.2f\n',totalTime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Now try manual, slower method
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
% Now let's place the floats along an isopycnal.
% isopycnalDeviation = wavemodel.ZetaAtTimePosition(0,x_float,y_float,z_float,interpolationMethod);
zIsopycnal3 = z_float;% + isopycnalDeviation;

% Iteratively place floats on the isopycnal surface. Overkill, probably.
for zLevel = 1:nLevels
    zLevelIndices = (zLevel-1)*N*N + (1:(N*N));
    for i=1:2
        rho = wavemodel.DensityAtTimePosition(0,x_float(zLevelIndices),y_float(zLevelIndices),zIsopycnal3(zLevelIndices),interpolationMethod);
        dRho = rho - mean(rho);
        dz = dRho * 9.81/(N0*N0*wavemodel.rho0);
        zIsopycnal3(zLevelIndices) = zIsopycnal3(zLevelIndices)+dz;
    end
    
    rho = wavemodel.DensityAtTimePosition(0,x_float(zLevelIndices),y_float(zLevelIndices),zIsopycnal3(zLevelIndices));
    dRho = rho - mean(rho);
    dz = dRho * 9.81/(N0*N0*wavemodel.rho0);
    fprintf('All floats are within %.2g meters of the isopycnal at z=%.1f meters\n',max(abs(dz)),z_float((zLevel-1)*N*N+1))
end

totalTime = toc;
fprintf('Total time %.2f\n',totalTime);