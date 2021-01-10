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

alpha = atan2(wm.L,wm.K);
g = wm.g;
omega = wm.Omega;
fOmega = wm.f0./omega;
Kh = wm.Kh;

ApU = (1/2)*(cos(alpha)+sqrt(-1)*fOmega.*sin(alpha));
ApV = (1/2)*(sin(alpha)-sqrt(-1)*fOmega.*cos(alpha));
ApN = -g*Kh./(2*omega);

AmU = (1/2)*(cos(alpha)-sqrt(-1)*fOmega.*sin(alpha));
AmV = (1/2)*(sin(alpha)+sqrt(-1)*fOmega.*cos(alpha));
AmN = g*Kh./(2*omega);

A0U = sqrt(-1)*wm.h.*(fOmega./omega) .* wm.L;
A0V = -sqrt(-1)*wm.h.*(fOmega./omega) .* wm.K; 
A0N = fOmega.^2;

Ap = ApU.*Ubar + ApV.*Vbar + ApN.*Nbar;
Am = AmU.*Ubar + AmV.*Vbar + AmN.*Nbar;
A0 = A0U.*Ubar + A0V.*Vbar + A0N.*Nbar;

Ap = sqrt(wm.h).*Ap;
Am = sqrt(wm.h).*Am;