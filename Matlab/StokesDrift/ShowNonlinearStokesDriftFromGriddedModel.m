%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjustable parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wavelength = 10e3; % wavelength in meters
j = 1; % vertical mode number
epsilon = 0.05; % nonlinearity parameter
maxOscillations = 5; % Total number of oscillations, in periods
stratification = 'constant';
z0 = [-10; -250; -625]; % initial particle positions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derived parameters & model setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rho, N2, zIn] = InternalModes.StratificationProfileWithName(stratification);
latitude = 0.0;
k = 2*pi/wavelength;
wavenumber = 4;
Nx = 128;
Ny = 2;
Nz = 513; % Must include end point to advect at the surface, so use 2^N + 1
Lx = wavenumber*wavelength;
Ly = Ny*Lx/Nx;
Lz = max(zIn)-min(zIn);
zDomain = linspace(min(zIn),max(zIn),Nz)';

wavemodel = InternalWaveModelArbitraryStratification([Lx, Ly, Lz], [Nx, Ny, Nz], rho, zDomain, j+1, latitude);

% Pull out the wave modes to sketch out the theory for later (and get the
% wave period).
wavemodel.internalModes.normalization = Normalization.uMax;
wavemodel.internalModes.nModes = floor(Nz/2);
[F_out,G_out,h_out,omega_out] = wavemodel.internalModes.ModesAtWavenumber(k);
F = F_out(:,j);
G = G_out(:,j);
h = h_out(j);
omega = omega_out(j);
U = epsilon*(omega/k);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the nonlinear correction
% We use much higher resolution modes than the model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zHR = linspace(min(zIn),max(zIn),5000)';
im = InternalModesWKBSpectral(rho,zIn,zHR,latitude);
im.normalization = Normalization.uMax;
[F_hr_out,G_hr_out,~,~] = im.ModesAtWavenumber(k);
F_hr = F_hr_out(:,j);
G_hr = G_hr_out(:,j);

f = -im.rho_zz .* G_hr .* G_hr / wavemodel.internalModes.rho0;

im.normalization = Normalization.kConstant;
[F_2k_out,G_2k_out,h_2k_out,~] = im.ModesAtWavenumber(2*k);
a = zeros(size(h_2k_out));
for i=1:size(G_2k_out,2)
    a(i) = trapz(zHR,f.*G_2k_out(:,i));
end
gamma = ((h*h_2k_out./(h-h_2k_out))).*a;

figure
subplot(1,2,1)
plot(F_hr,zHR), hold on
plot(F_2k_out*((h./h_2k_out).*gamma).',zHR)
title('U-mode')
subplot(1,2,2)
plot(G_hr,zHR), hold on
plot(G_2k_out*gamma.',zHR)
plot(h* (im.rho_zz./im.rho_z) .* G_hr .* G_hr, zHR)
title('W-mode')
legend('O(\epsilon)','O(\epsilon^2)','O(\epsilon^2)')

wavemodel.internalModes.normalization = Normalization.kConstant;
[F_2k_out,G_2k_out,~,~] = wavemodel.internalModes.ModesAtWavenumber(2*k);

X = wavemodel.X;
Y = wavemodel.Y;
Z = wavemodel.Z;
RhoBarDz = -(wavemodel.rho0/9.81)*wavemodel.N2AtDepth(Z);

Phi = repmat(permute((F_2k_out*(((h./h_2k_out).*gamma).')),[3 2 1]),size(X,1),size(X,2));
Gamma = repmat(permute((G_2k_out*(gamma.')),[3 2 1]),size(X,1),size(X,2));
Rho_zzTerm = repmat(permute((h*wavemodel.internalModes.rho_zz.*G.*G),[3 2 1]),size(X,1),size(X,2));

u_correction = (U * epsilon / 4) * cos(2*k*X) .* Phi;
w_correction = (U*k*h*epsilon/2) * sin(2*k*X) .* Gamma;
rho_correction = (h*epsilon*epsilon/4) * cos(2*k*X) .* ( Rho_zzTerm +  RhoBarDz .* Gamma );

period = wavemodel.InitializeWithPlaneWave(wavenumber,0,j,U,1);

[u_i,v_i] = wavemodel.VelocityFieldAtTime(0);
w_i = wavemodel.VerticalFieldsAtTime(0);
rho_i = wavemodel.DensityAtTime(0);

figure, plot(squeeze(u_i(1,1,:)),wavemodel.z), hold on
plot(squeeze(u_correction(1,1,:)),wavemodel.z)

figure, plot(squeeze(w_i(Nx/2,1,:)),wavemodel.z), hold on
plot(squeeze(w_correction(Nx/2,1,:)),wavemodel.z)

figure, plot(squeeze(rho_i(1,1,:))-wavemodel.RhoBarAtDepth(wavemodel.z),wavemodel.z), hold on
plot(squeeze(rho_correction(1,1,:)),wavemodel.z)

return

t = (0:(1/40):maxOscillations)'*period;
fprintf('The wave period is set to %.1f hours.\n', period/3600)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Velocity field & numerical time stepping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = @(t, x) wavemodel.VelocityAtTimePositionVector(t,x);
x = zeros(length(t),length(z0));
y = zeros(length(t),length(z0));
z = zeros(length(t),length(z0));
b = zeros(length(t),length(z0));

for i=1:length(z0)
    % Using ode113 for extremely high error tolerances
    options = odeset('RelTol',1e-5,'AbsTol',1e-5); % overkill
    [~, X] = ode113(u,t,[0 0 z0(i)],options);
    
    x(:,i) = X(:,1);
    y(:,i) = X(:,2);
    z(:,i) = X(:,3);
%     b(:,i) = rho(t,[x(:,i),z(:,i)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Theoretical Stokes drift
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zTheory = linspace(-1000,0,1000)';
U1 = (U*U*k)/(2*omega)*(interp1(zDomain,F,zTheory)).^2;
U2 = -(U*U*k)/(2*omega)*(h/g*(N2(zTheory)-omega*omega)).*(interp1(zDomain,G,zTheory)).^2;

figure
plot(-(U1+U2)*max(t),zTheory)
hold on
plot(x,z,'LineWidth',2)

