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

wvm = WaveVortexModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);

if shouldUseGMSpectrum == 1
    wvm.initWithGMSpectrum(1.0);
    maxTime = 60*60; %2*pi/wavemodel.f;
    period = 2*pi/wvm.N0;
    [u,v] = wvm.VelocityFieldAtTime(0.0);
    U = max(max(max( sqrt(u.*u + v.*v) )));
else
    j0 = 1; % j=1..nModes, where 1 indicates the 1st baroclinic mode
    U = 0.01; % m/s
    sign = 1;
    phi = 0;
    k0 = 1;
    l0 = 0;
    alpha = atan2(l0,k0);
    k = 2*pi*sqrt(k0^2 + l0^2)/Lx;
    
    period = wvm.initWithWaveModes(k0,l0,j0,phi,U,sign);
    maxTime = period;
end


% initial positions, center of the horizontal domain doesn't deal with boundaries
p0 = [Lx/2, Ly/2, -1*Lz/4];

% iniital position near horizontal boundary, requires dealing with the periodicity of the boundary
% p0 = [0, 0, -0*Lz/4];

f = @(t,y) wvm.VelocityAtTimePositionVector(t,y,'spline');

% Let's do fixed step size integrator.
cfl = 0.25;
advectiveDT = cfl*(wvm.x(2)-wvm.x(1))/U;
oscillatoryDT = period/8;
if advectiveDT < oscillatoryDT
    fprintf('Using the advective dt: %.2f\n',advectiveDT);
    deltaT = advectiveDT;
else
    fprintf('Using the oscillatory dt: %.2f\n',oscillatoryDT);
    deltaT = oscillatoryDT;
end

t_in = (0:deltaT:maxTime)'; %(0:60*ceil(deltaT/60):maxTime)';
if t_in(end) < period
    t_in(end+1) = period;
end

% Try ode2, ode3, & ode4, depending on required accuracy.
tic
p2 = ode2(f,t_in, p0');
t2 = toc;
x2 = p2(:,1);
y2 = p2(:,2);
z2 = p2(:,3);

tic
p3 = ode3(f,t_in, p0');
t3 = toc;
x3 = p3(:,1);
y3 = p3(:,2);
z3 = p3(:,3);

tic
p4 = ode4(f,t_in, p0');
t4 = toc;
x4 = p4(:,1);
y4 = p4(:,2);
z4 = p4(:,3);

tic
integrator = Integrator(f,p0',deltaT);
pI = integrator.IntegrateAlongDimension(t_in);
tI = toc;
xI = pI(:,1);
yI = pI(:,2);
zI = pI(:,3);

% profile on
tic
[~,p45] = ode45(f,t_in, p0); % ,odeset('RelTol',1e-11,'AbsTol',1e-11)
t45 = toc;
x45 = p45(:,1);
y45 = p45(:,2);
z45 = p45(:,3);
% profile viewer

% Second, let's do the adaptive time-stepping integrator
if shouldUseGMSpectrum == 1
    [omega, alpha, mode, phi, A] = wvm.waveModesFromWaveCoefficients();
    wvm.initWithWaveModes(1,1,1,0,0.0,1);
    wvm.setExternalWavesWithFrequencies(omega, alpha, mode, phi, A,Normalization.kConstant);
else
    wvm.initWithWaveModes(1,1,1,0,0.0,1);
    k0 = k*cos(alpha);
    l0 = k*sin(alpha);
    omega = wvm.setExternalWavesWithWavenumbers(k0,l0,j0,phi,U,Normalization.uMax);
end

% This has to be repeated, to capture the new wavemodel reference.
% f = @(t,y) wavemodel.VelocityAtTimePositionVector(t,y);

f = @(t,y) wvm.VelocityAtTimePositionVector(t,y,'exact');

tic
[t,p] = ode45(f,t_in, p0,odeset('RelTol',1e-11,'AbsTol',1e-8)); % ,odeset('RelTol',1e-11,'AbsTol',1e-11)
tAdaptive = toc;
x = p(:,1);
y = p(:,2);
z = p(:,3);
stokes_x = x(end)-x(1);
fprintf('Adaptive time step took t=%.2f seconds.\n',tAdaptive);

figure
plot(x,y,'LineWidth',2), hold on, plot(x45,y45), plot(x2,y2), plot(x3,y3), plot(x4,y4)
legend('adaptive time step (exact)', 'adaptive time step (gridded)', 'fixed time step (ode2)', 'fixed time step (ode3)', 'fixed time step (ode4)')

% This is a measure of relative error useful for position.
error = @(u,u_unit) max( [max(max(max(abs(u-u_unit)/(max(max(max( u_unit ))) - min(min(min( u_unit )))) ))), 1e-15]);
kappa_h_effective = @(x,x_unit,y,y_unit) 0.25*((x(end)-x_unit(end)).^2 + (y(end)-y_unit(end)).^2)/(t(end));
% kappa_h_effective = @(x,x_unit,y,y_unit) max(0.25*((x(2:end)-x_unit(2:end)).^2 + (y(2:end)-y_unit(2:end)).^2)./(t(2:end)));
kappa_z_effective = @(z,z_unit) 0.5*(z(end)-z_unit(end)).^2/(t(end));

x_error = error(x45,x);
y_error = error(y45,y);
z_error = error(z45,z);
kappa_h_eff45 = kappa_h_effective(x45,x,y45,y);
kappa_z_eff45 = kappa_z_effective(z45,z);
stokes_x45 = x45(end)-x45(1);
fprintf('The ode45 solution for (x,y,z) matches the exact RK45 solution to 1 part in (10^%d, 10^%d, 10^%d). K_eff = (%.2g, %.2g). Stokes ratio: %.2f. t=%.2f\n', round((log10(x_error))), round((log10(y_error))), round((log10(z_error))), kappa_h_eff45, kappa_z_eff45,stokes_x45/stokes_x, t45);

x_error = error(x2,x);
y_error = error(y2,y);
z_error = error(z2,z);
kappa_h_eff2 = kappa_h_effective(x2,x,y2,y);
kappa_z_eff2 = kappa_z_effective(z2,z);
stokes_x2 = x2(end)-x2(1);
fprintf('The ode2 solution for (x,y,z) matches the RK45 solution to 1 part in (10^%d, 10^%d, 10^%d). K_eff = (%.2g, %.2g). Stokes ratio: %.2f. t=%.2f\n', round((log10(x_error))), round((log10(y_error))), round((log10(z_error))), kappa_h_eff2, kappa_z_eff2,stokes_x2/stokes_x, t2);

x_error = error(x3,x);
y_error = error(y3,y);
z_error = error(z3,z);
kappa_h_eff3 = kappa_h_effective(x3,x,y3,y);
kappa_z_eff3 = kappa_z_effective(z3,z);
stokes_x3 = x3(end)-x3(1);
fprintf('The ode3 solution for (x,y,z) matches the RK45 solution to 1 part in (10^%d, 10^%d, 10^%d). K_eff = (%.2g, %.2g). Stokes ratio: %.2f. t=%.2f\n', round((log10(x_error))), round((log10(y_error))), round((log10(z_error))), kappa_h_eff3, kappa_z_eff3,stokes_x3/stokes_x, t3);

x_error = error(x4,x);
y_error = error(y4,y);
z_error = error(z4,z);
kappa_z_eff4 = kappa_z_effective(z4,z);
kappa_h_eff4 = kappa_h_effective(x4,x,y4,y);
stokes_x4 = x4(end)-x4(1);
fprintf('The ode4 solution for (x,y,z) matches the RK45 solution to 1 part in (10^%d, 10^%d, 10^%d). K_eff = (%.2g, %.2g). Stokes ratio: %.2f. t=%.2f\n', round((log10(x_error))), round((log10(y_error))), round((log10(z_error))), kappa_h_eff4, kappa_z_eff4,stokes_x4/stokes_x, t4);

