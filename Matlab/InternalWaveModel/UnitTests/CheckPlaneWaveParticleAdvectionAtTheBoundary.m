Lx = 15e3;
Ly = 15e3;
Lz = 5000;

Nx = 64;
Ny = 64;
Nz = 65; % Must include end point to advect at the surface, so use 2^N + 1

latitude = 31;
N0 = 5.2e-3; % Choose your stratification

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shouldUseGMSpectrum = 0;

wavemodel = InternalWaveModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);

if shouldUseGMSpectrum == 1
    wavemodel.InitializeWithGMSpectrum(1.0);
    maxTime = 60*60; %2*pi/wavemodel.f0;
    period = 2*pi/wavemodel.N0;
    [u,v] = wavemodel.VelocityFieldAtTime(0.0);
    U = max(max(max( sqrt(u.*u + v.*v) )));
else
    j0 = 1; % j=1..nModes, where 1 indicates the 1st baroclinic mode
    U = 0.01; % m/s
    sign = 1;
    phi = 0;
    k0 = 16;
    l0 = 0;
    alpha = atan2(l0,k0);
    k = 2*pi*sqrt(k0^2 + l0^2)/Lx;
    
    period = wavemodel.InitializeWithPlaneWave(k0,l0,j0,U,sign);
    maxTime = period;
end

deltaT=maxTime/21;
t_in = (0:deltaT:maxTime)';

dx = wavemodel.x(2) - wavemodel.x(1);
dy = wavemodel.y(2) - wavemodel.y(1);

% initial positions, center of the horizontal domain doesn't deal with boundaries
p0 = [Lx/2, Ly/2, -1*Lz/4];

% iniital position near horizontal boundary, requires dealing with the periodicity of the boundary
% p0 = [0, 0, -1*Lz/4];

f = @(t,y) wavemodel.VelocityAtTimePositionVector(t,y,'exact');
f_spline = @(t,y) wavemodel.VelocityAtTimePositionVector(t,y,'spline');

[t,p] = ode45(f,t_in, p0);
x = p(:,1);
y = p(:,2);
z = p(:,3);

[t,p] = ode45(f_spline,t_in, p0);
xSpline = p(:,1);
ySpline = p(:,2);
zSpline = p(:,3);

figure
plot(x,y), hold on, plot(xSpline,ySpline)

% This is a measure of relative error useful for position.
error = @(u,u_unit) max( [max(max(max(abs(u-u_unit)/(max(max(max( u_unit ))) - min(min(min( u_unit )))) ))), 1e-15]);

x_error = error(xSpline,x);
y_error = error(ySpline,y);
z_error = error(zSpline,z);

fprintf('The ode45 solution for (x,y,z) matches the exact RK45 solution to 1 part in (10^%d, 10^%d, 10^%d).\n', round((log10(x_error))), round((log10(y_error))), round((log10(z_error))) );
