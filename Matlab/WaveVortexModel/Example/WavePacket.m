%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 256;
aspectRatio = 1/2;

Lx = 100e3;
Ly = aspectRatio*Lx;
Lz = 1300;

Nx = N;
Ny = aspectRatio*N;
Nz = N+1; % 2^n + 1 grid points, to match the Winters model, but 2^n ok too.

latitude = 31;
N0 = 5.2e-3/2; % Choose your stratification 7.6001e-04

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wvm = WaveVortexModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);

cg_x = wvm.cg_x;
cg_y = wvm.cg_y;
cg_z = wvm.cg_z;

a = squeeze(cg_z(:,10,:));
figure, pcolor(a), shading interp, colorbar('eastoutside'), xlog
% b = squeeze(cg_x(:,1,:));
% figure, pcolor(b), shading interp, colorbar('eastoutside'), xlog

[X,Y,Z] = ndgrid(wvm.x,wvm.y,wvm.z);
[K,L,J] = ndgrid(wvm.k,wvm.l,wvm.j);
alpha = atan2(L,K);
K2 = K.*K + L.*L;
Kh = sqrt(K2);

Lh = Lx/16;
Lv = Lz/8;
x0 = Lx/2;
y0 = Ly/2;
z0 = -Lz/2;
eta0 = 100*exp( -((X-x0).^2 + (Y-y0).^2)/(Lh)^2  - ((Z-z0).^2)/(Lv)^2 ).*sin(X/(Lh/32)+Z/(Lv/16));

eta0_bar = wvm.transformFromSpatialDomainWithG(eta0);
A_plus = eta0_bar ./ wvm.NAp;
A_plus(isnan(A_plus)) = 0;
A_plus(isinf(A_plus)) = 0;
% A_plus = wvm.ApN .* eta0_bar;
A_plus(K < 0) = 0;
wvm.Ap = A_plus;

[u, v, w, rho_prime, eta, p_wave]= wvm.VariableFieldsAtTime(75*3600, 'u', 'v', 'w', 'rho_prime', 'eta', 'p');

maxU = max(max(max(abs(u))));
maxV = max(max(max(abs(v))));
maxW = max(max(max(abs(w))));
fprintf('Maximum fluid velocity (u,v,w)=(%.2f,%.2f,%.2f) cm/s\n',100*maxU,100*maxV,100*maxW);

cutoff = max(abs(A_plus(:)))/20;
figure
subplot(1,2,1)
histogram( 100*cg_x(abs(A_plus(:))>cutoff), 'Normalization', 'cdf' )
subplot(1,2,2)
histogram( 100*cg_z(abs(A_plus(:))>cutoff), 'Normalization', 'cdf' )

dispvar = eta;
figure
subplot(2,1,1)
pcolor(wvm.x/1000,wvm.z,squeeze(dispvar(:,Ny/2,:))'),shading flat
subplot(2,1,2)
pcolor(wvm.x/1000,wvm.y/1000,squeeze(dispvar(:,:,floor(Nz/2)))'),shading flat, axis equal

[u, v, w, rho_prime, eta, p_wave]= wvm.VariableFieldsAtTime(250*3600, 'u', 'v', 'w', 'rho_prime', 'eta', 'p');
dispvar = eta;
figure
subplot(2,1,1)
pcolor(wvm.x/1000,wvm.z,squeeze(dispvar(:,Ny/2,:))'),shading interp
subplot(2,1,2)
pcolor(wvm.x/1000,wvm.y/1000,squeeze(dispvar(:,:,floor(Nz/2)))'),shading interp, axis equal
return;

U = .2;
% boussinesq.setWaveModes(0,0,1,U,1); 
[omega,k,l] = wvm.initWithWaveModes(10, 0, 1, 0, U, 1)
period = 2*pi/omega; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the integrator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = period/50;
nT=5*50;
nTrajectories = 101;
totalEnergy = zeros(nT,1);
totalSpectralEnergy = zeros(nT,1);
totalEnergy(1) = wvm.totalEnergy;
totalSpectralEnergy(1) = wvm.totalSpectralEnergy;
x = zeros(nT,nTrajectories); y = zeros(nT,nTrajectories); z = zeros(nT,nTrajectories);
x(1,:) = Lx/2*ones(1,nTrajectories); y(1,:) = Ly/2*ones(1,nTrajectories); z(1,:) = linspace(-Lz,0,nTrajectories);

integrator = ArrayIntegrator(@(t,y0) wvm.NonlinearFluxWithParticlesAtTimeArray(t,y0),{wvm.Ap,wvm.Am,wvm.A0,x(1,:),y(1,:),z(1,:)},dt);

% profile on
for i=2:nT
%    integrator.currentY = wvm.Y;
   integrator.IncrementForward();
   wvm.Ap = integrator.currentY{1};
   wvm.Am = integrator.currentY{2};
   wvm.A0 = integrator.currentY{3};
   totalEnergy(i) = wvm.totalEnergy;
   totalSpectralEnergy(i) = wvm.totalSpectralEnergy;
   x(i,:) = integrator.currentY{4};
   y(i,:) = integrator.currentY{5};
   z(i,:) = integrator.currentY{6};
%    if mod(i,10)==0
       wvm.summarizeEnergyContent();
%    end
end
% profile viewer