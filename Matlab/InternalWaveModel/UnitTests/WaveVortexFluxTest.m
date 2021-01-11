%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CheckWavenumberCutoff
%
% Tests the API for limiting the the horizontal wavenumber used to
% initialize waves.
%
% Jeffrey J. Early
% jeffrey@jeffreyearly.com
%
% February 8, 2018      Version 1.0


N = 16;
aspectRatio = 1;

Lx = 100e3;
Ly = aspectRatio*100e3;
Lz = 5000;

Nx = N;
Ny = aspectRatio*N;
Nz = N+1; % 2^n + 1 grid points, to match the Winters model, but 2^n ok too.

latitude = 31;
N0 = 5.2e-3; % Choose your stratification 7.6001e-04

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wm = InternalWaveModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);
wm.InitializeWithGMSpectrum(1.0)
[u,v,w,eta]=wm.VariableFieldsAtTime(0,'u','v','w','zeta');

Ubar = wm.TransformFromSpatialDomainWithF(u);
Vbar = wm.TransformFromSpatialDomainWithF(v);
Nbar = wm.TransformFromSpatialDomainWithG(eta);

K = wm.K;
L = wm.L;
alpha = atan2(wm.L,wm.K);
g = wm.g;
omega = wm.Omega;
fOmega = wm.f0./omega;
Kh = wm.Kh;

ApU = wm.MakeHermitian((1/2)*(cos(alpha)+sqrt(-1)*fOmega.*sin(alpha)));
ApV = wm.MakeHermitian((1/2)*(sin(alpha)-sqrt(-1)*fOmega.*cos(alpha)));
ApN = wm.MakeHermitian(-g*Kh./(2*omega));

AmU = wm.MakeHermitian((1/2)*(cos(alpha)-sqrt(-1)*fOmega.*sin(alpha)));
AmV = wm.MakeHermitian((1/2)*(sin(alpha)+sqrt(-1)*fOmega.*cos(alpha)));
AmN = wm.MakeHermitian(g*Kh./(2*omega));

A0U = wm.MakeHermitian(sqrt(-1)*wm.h.*(fOmega./omega) .* wm.L);
A0V = wm.MakeHermitian(-sqrt(-1)*wm.h.*(fOmega./omega) .* wm.K); 
A0N = wm.MakeHermitian(fOmega.^2);

Ap = ApU.*Ubar + ApV.*Vbar + ApN.*Nbar;
Am = AmU.*Ubar + AmV.*Vbar + AmN.*Nbar;
A0 = A0U.*Ubar + A0V.*Vbar + A0N.*Nbar;

bm = Boussinesq3DConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);
[App,Amm,A00] = bm.Project(u,v,eta);

% Ap = sqrt(wm.h).*Ap;
% Am = sqrt(wm.h).*Am;