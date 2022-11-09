N = 64;
aspectRatio = 1;

Lx = 100e3;
Ly = aspectRatio*100e3;
Lz = 1300;

Nx = N;
Ny = aspectRatio*N;
Nz = N+1; % 2^n + 1 grid points, to match the Winters model, but 2^n ok too.

latitude = 31;
N0 = 5.2e-3; % Choose your stratification 7.6001e-04
rho0 = 1025;

outputfile = '/Users/jearly/Data/InternalWaveSimulation/igw-simulation.nc';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wvt = WVTransformConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], N0,latitude=latitude);
rng(1)
% wvt.initWithGMSpectrum(1.0)
wvt.initWithRandomFlow();
wvt.removeEnergyFromAliasedModes();
wvt.nonlinearFluxOperation = BoussinesqConstantN(wvt,shouldAntialias=1);

addpath('/Users/jearly/Documents/ProjectRepositories/GLOceanKit-copy/Matlab/WaveVortexModel');
wvm = WaveVortexModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0, rho0, 'shouldAntiAlias', 1);
rng(1)
% wvm.InitializeWithGMSpectrum(1.0);
wvm.Ap = wvt.Ap;
wvm.Am = wvt.Am;
wvm.A0 = wvt.A0;
% wvm.clearEnergyFromAliasedModes;

error = @(u,u_unit) max( [max(max(max(abs(u-u_unit)/max( [max(max(max( abs(u) ))), 1e-15] )))), 1e-15]);
error2 = @(u,u_unit) abs((u-u_unit))./(max(max(max(abs(u_unit)))));
error3 = @(u,u_unit) abs((u-u_unit)./max(abs(u_unit),1e-15));

Ap_error = error3(wvt.Ap,wvm.Ap);
Am_error = error3(wvt.Am,wvm.Am);
A0_error = error3(wvt.A0,wvm.A0);
fprintf('\tAp error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(Ap_error)))))));
fprintf('\tAm error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(Am_error)))))));
fprintf('\tA0 error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(A0_error)))))));

[wvtFp,wvtFm,wvtF0] = wvt.nonlinearFlux;
[wvmFp,wvmFm,wvmF0] = wvm.NonlinearFluxAtTime(0,wvm.Ap,wvm.Am,wvm.A0);

Fp_error = error3(wvtFp,wvmFp);
Fm_error = error3(wvtFm,wvmFm);
F0_error = error3(wvtF0,wvmF0);

fprintf('\tFp error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(Fp_error)))))));
fprintf('\tFm error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(Fm_error)))))));
fprintf('\tF0 error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(F0_error)))))));

dFp = wvtFp-wvmFp;
dFm = wvtFm-wvmFm;
dF0 = wvtF0-wvmF0;

[wvtEp,wvtEm,wvtE0] = wvt.energyFlux;
[wvmEp,wvmEm,wvmE0] = wvm.EnergyFluxAtTime(0,wvm.Ap,wvm.Am,wvm.A0);

Ep_error = error3(wvtEp,wvmEp);
Em_error = error3(wvtEm,wvmEm);
E0_error = error3(wvtE0,wvmE0);

fprintf('\tEp error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(Ep_error)))))));
fprintf('\tEm error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(Em_error)))))));
fprintf('\tE0 error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(E0_error)))))));

dEp = wvtEp-wvmEp;
dEm = wvtEm-wvmEm;
dE0 = wvtE0-wvmE0;

% Nonlinear terms with new wvt.
Ubar = wvt.UAp.*wvt.Ap + wvt.UAm.*wvt.Am + wvt.UA0.*wvt.A0;
Vbar = wvt.VAp.*wvt.Ap + wvt.VAm.*wvt.Am + wvt.VA0.*wvt.A0;
Wbar = wvt.WAp.*wvt.Ap + wvt.WAm.*wvt.Am;
Nbar = wvt.NAp.*wvt.Ap + wvt.NAm.*wvt.Am + wvt.NA0.*wvt.A0;
[U1,Ux1,Uy1,Uz1] = wvt.transformToSpatialDomainWithFAllDerivatives(Ubar);
[V1,Vx1,Vy1,Vz1] = wvt.transformToSpatialDomainWithFAllDerivatives(Vbar);
W1 = wvt.transformToSpatialDomainWithG(Wbar);
[~,ETAx1,ETAy1,ETAz1] = wvt.transformToSpatialDomainWithGAllDerivatives(Nbar);
uNL1 = -U1.*Ux1 - V1.*Uy1 - W1.*Uz1;
vNL1 = -U1.*Vx1 - V1.*Vy1 - W1.*Vz1;
nNL1 = -U1.*ETAx1 - V1.*ETAy1 - W1.*ETAz1;

% Ubar = wvm.UAp.*wvt.Ap + wvm.UAm.*wvt.Am + wvm.UA0.*wvt.A0;
% Vbar = wvm.VAp.*wvt.Ap + wvm.VAm.*wvt.Am + wvm.VA0.*wvt.A0;
% Wbar = wvm.WAp.*wvt.Ap + wvm.WAm.*wvt.Am;
% Nbar = wvm.NAp.*wvt.Ap + wvm.NAm.*wvt.Am + wvm.NA0.*wvt.A0;
% [U,Ux,Uy,Uz] = wvm.TransformToSpatialDomainWithFAllDerivatives(Ubar);
% [V,Vx,Vy,Vz] = wvm.TransformToSpatialDomainWithFAllDerivatives(Vbar);
% W = wvm.TransformToSpatialDomainWithG(Wbar);
% [~,ETAx,ETAy,ETAz] = wvm.TransformToSpatialDomainWithGAllDerivatives(Nbar);
% uNL1 = -U.*Ux - V.*Uy - W.*Uz;
% vNL1 = -U.*Vx - V.*Vy - W.*Vz;
% nNL1 = -U.*ETAx - V.*ETAy - W.*ETAz;

% Nonlinear terms with old wvm.
Ubar = wvm.UAp.*wvm.Ap + wvm.UAm.*wvm.Am + wvm.UA0.*wvm.A0;
Vbar = wvm.VAp.*wvm.Ap + wvm.VAm.*wvm.Am + wvm.VA0.*wvm.A0;
Wbar = wvm.WAp.*wvm.Ap + wvm.WAm.*wvm.Am;
Nbar = wvm.NAp.*wvm.Ap + wvm.NAm.*wvm.Am + wvm.NA0.*wvm.A0;
[U,Ux,Uy,Uz] = wvm.TransformToSpatialDomainWithFAllDerivatives(Ubar);
[V,Vx,Vy,Vz] = wvm.TransformToSpatialDomainWithFAllDerivatives(Vbar);
W = wvm.TransformToSpatialDomainWithG(Wbar);
[~,ETAx,ETAy,ETAz] = wvm.TransformToSpatialDomainWithGAllDerivatives(Nbar);
uNL2 = -U.*Ux - V.*Uy - W.*Uz;
vNL2 = -U.*Vx - V.*Vy - W.*Vz;
nNL2 = -U.*ETAx - V.*ETAy - W.*ETAz;

u_error = error3(U,U1);
v_error = error3(V,V1);
w_error = error3(W,W1);
ux_error = error3(Ux,Ux1);
vx_error = error3(Vx,Vx1);
uy_error = error3(Uy,Uy1);
vy_error = error3(Vy,Vy1);
uz_error = error3(Uz,Uz1);
vz_error = error3(Vz,Vz1);
fprintf('\tu error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(u_error)))))));
fprintf('\tv error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(v_error)))))));
fprintf('\tw error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(w_error)))))));
fprintf('\tux error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(ux_error)))))));
fprintf('\tvx error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(vx_error)))))));
fprintf('\tuy error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(uy_error)))))));
fprintf('\tvy error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(vy_error)))))));
fprintf('\tuz error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(uz_error)))))));
fprintf('\tvz error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(vz_error)))))));

uNL_error = error3(uNL1,uNL2);
vNL_error = error3(vNL1,vNL2);
nNL_error = error3(nNL1,nNL2);

fprintf('\tuNL error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(uNL_error)))))));
fprintf('\tvNL error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(vNL_error)))))));
fprintf('\tnNL error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(nNL_error)))))));

% Now apply the operator S^{-1} and then T_\omega^{-1}
uNLbar = wvt.transformFromSpatialDomainWithF(uNL1);
vNLbar = wvt.transformFromSpatialDomainWithF(vNL1);
nNLbar = wvt.transformFromSpatialDomainWithG(nNL1);

Fp1 = wvm.EMAp .* (wvt.ApU.*uNLbar + wvt.ApV.*vNLbar + wvt.ApN.*nNLbar);
Fm1 = wvm.EMAm .* (wvt.AmU.*uNLbar + wvt.AmV.*vNLbar + wvt.AmN.*nNLbar);
F01 = wvm.EMA0 .* (wvt.A0U.*uNLbar + wvt.A0V.*vNLbar + wvt.A0N.*nNLbar);

uNLbar = wvm.TransformFromSpatialDomainWithF(uNL2);
vNLbar = wvm.TransformFromSpatialDomainWithF(vNL2);
nNLbar = wvm.TransformFromSpatialDomainWithG(nNL2);

Fp2 = wvm.EMAp .* (wvm.ApU.*uNLbar + wvm.ApV.*vNLbar + wvm.ApN.*nNLbar);
Fm2 = wvm.EMAm .* (wvm.AmU.*uNLbar + wvm.AmV.*vNLbar + wvm.AmN.*nNLbar);
F02 = wvm.EMA0 .* (wvm.A0U.*uNLbar + wvm.A0V.*vNLbar + wvm.A0N.*nNLbar);

Fp2_error = error3(Fp1,Fp2);
Fm2_error = error3(Fm1,Fm2);
F02_error = error3(F01,F02);

fprintf('\tFp error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(Fp2_error)))))));
fprintf('\tFm error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(Fm2_error)))))));
fprintf('\tF0 error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(F02_error)))))));
