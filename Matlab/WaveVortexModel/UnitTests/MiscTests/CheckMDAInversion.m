N0 = 3*2*pi/3600; % buoyancy frequency at the surface, radians/seconds
L_gm = 1300; % thermocline exponential scale, meters
Lz = 4000; % depth of the ocean
N2 = @(z) N0*N0*exp(2*z/L_gm);
% N2 = @(z) N0*N0*ones(size(z));

z = linspace(-Lz,0,1000)';
nModes = 64;
nPoints = nModes+1;
im = InternalModesSpectral(N2=N2,zIn=[-Lz 0],zOut=z,latitude=33,nModes=nModes);
im.normalization = Normalization.wMax;


%% Quadrature grid for geostrophic modes
im.upperBoundary = UpperBoundary.rigidLid;
im.lowerBoundary = LowerBoundary.freeSlip;
z_g = im.GaussQuadraturePointsForModesAtFrequency(nPoints,0);
im_g = InternalModesSpectral(N2=N2,zIn=[-Lz 0],zOut=z_g,latitude=33,nModes=nModes);
im_g.upperBoundary = UpperBoundary.rigidLid;
im_g.lowerBoundary = LowerBoundary.freeSlip;
im_g.normalization = Normalization.kConstant;
[Fg,Gg,hg] = im_g.ModesAtFrequency(0);

Pg = max(abs(Fg),[],1);
FgNorm = Fg./Pg;

%% Quadrature grid for mda modes
z_mda = im.GaussQuadraturePointsForMDAModes(nPoints);
im_mda = InternalModesSpectral(N2=N2,zIn=[-Lz 0],zOut=z_mda,latitude=33,nModes=nModes);
im_mda.normalization = Normalization.kConstant;
[Fmda,Gmda,hmda] = im_mda.MDAModes();

Qg = max(abs(Gmda),[],1);
GmdaNorm = Gmda./Qg;


%% MDA modes on the geostrophic quadrature grid
im_mda_g = InternalModesSpectral(N2=N2,zIn=[-Lz 0],zOut=z_g,latitude=33,nModes=nModes+1);
im_mda_g.normalization = Normalization.geostrophic;
[FmdaInv,GmdaInv,hmda_g] = im_mda_g.MDAModes();

Qg = max(abs(GmdaInv),[],1);
GmdaNorm = GmdaInv./Qg;
cond(GmdaNorm(2:end,:))

Pg = max(abs(FmdaInv),[],1);
FmdaNorm = FmdaInv./Pg;
cond(FmdaNorm(1:end-1,:))

A = FmdaInv;
A(end,:) = GmdaInv(end,:);
cond(A)
Ainv = inv(A);

% So what's the problem?
% We make some nice MDA modes with G, project onto those, and because the
% G(z) = 1 mode exists to make things complete, these modes cannot keep
% G(0)=0 without it. If I remove the mode from the set, then how do I
% invert?
% So, okay, you just need modes with 

%%
wvt = WVTransformHydrostatic([800e3, 800e3, 4000],[64, 64, nPoints], N2=N2,latitude=30);
U = -0.30; % m/s
Le = 120e3;
He = wvt.Lz/5;
x0 = (1/2)*max(wvt.x); y0=max(wvt.y)/2;
% Finv = wvt.FinvMatrix;
% Ginv = wvt.GinvMatrix;
% F2 = Finv(:,2); G2 = Ginv(:,2); h2 = wvt.h(2);
% w = @(z) (wvt.Lz/2 - abs(wvt.Lz/2 + z)).*interp1(wvt.z,F2,z) + h2*interp1(wvt.z,G2,z);
% A = min(w(wvt.z));
% w = @(z) ((wvt.Lz/2 - abs(wvt.Lz/2 + z)).*interp1(wvt.z,F2,z) + h2*interp1(wvt.z,G2,z))/A;

%%
pbar = @(z) (pi*Le*Le/(wvt.Lx*wvt.Ly))*U*(Le/sqrt(2))*exp(1/2)*exp(-(z/He).^2 );
psi = @(x,y,z) U*(Le/sqrt(2))*exp(1/2)*exp(-((x-x0)/Le).^2 -((y-y0)/Le).^2 -(z/He).^2 ) - pbar(z);
% psi = @(x,y,z) U*(Le/sqrt(2))*exp(1/2)*exp(-((x-x0)/Le).^2 -((y-y0)/Le).^2 -(z/He).^2 ) - (pi*Le*Le/(wvt.Lx*wvt.Ly))*U*(Le/sqrt(2))*exp(1/2)*exp(-(z/He).^2 );

[X,Y,Z] = wvt.xyzGrid;
psi_z = wvt.diffZF((wvt.f/wvt.g)*psi(X,Y,Z));
psibar = squeeze(mean(mean(psi(X,Y,Z),1),2));
psizbar = squeeze(mean(mean(psi_z,1),2));

figure, tiledlayout(1,2), nexttile, plot(psibar,wvt.z), hold on, plot(pbar(wvt.z),wvt.z), nexttile, plot(psizbar,wvt.z)
figure, plot(psibar-pbar(wvt.z),wvt.z)

%%
figure
plot

return

%%
% If I have an MDA mode on a geostrophic quadrature grid, I can evaluate
% which geostrophic modes each MDA mode projects onto. How do I evaluate
% whether or not the MDA modes are being resolved?

% Turns out, the grids are not that far apart.
figure, plot(cat(2,zeros(size(z_g)),ones(size(z_g))).',[z_g,z_mda].')
figure, scatter(diff(z_g),diff(z_mda))
%%
modes=1:4;
figure
tl = tiledlayout(2,2,"TileSpacing","compact");

nexttile
ax = gca; ax.ColorOrderIndex = 1;
plot(GmdaInv(:,modes),im_g.z, 'LineWidth', 2), hold on
plot([0 0],[-Lz 0],LineWidth=1,Color=0*[1 1 1])
nexttile
ax = gca; ax.ColorOrderIndex = 1;
plot(FmdaInv(:,modes),im_g.z, 'LineWidth', 2), hold on
plot([0 0],[-Lz 0],LineWidth=1,Color=0*[1 1 1])
nexttile
ax = gca; ax.ColorOrderIndex = 1;
plot(Gmda(:,modes),im_mda.z, 'LineWidth', 2), hold on
plot([0 0],[-Lz 0],LineWidth=1,Color=0*[1 1 1])
nexttile
ax = gca; ax.ColorOrderIndex = 1;
plot(Fmda(:,modes),im_mda.z, 'LineWidth', 2), hold on
plot([0 0],[-Lz 0],LineWidth=1,Color=0*[1 1 1])