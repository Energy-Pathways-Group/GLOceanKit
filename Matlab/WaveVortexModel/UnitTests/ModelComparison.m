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

% The issue is clearly with the transformation from (Ap,Am,A0) ->
% (uNL,vNL,nNL) using the new transform. But what aspect of that
% transformation?

[Ubar1,Vbar1,Wbar1,Nbar1] = physicalVarsFromWaveVortexVars(wvm,wvt.Ap,wvt.Am,wvt.A0);
[uNL1,vNL1,nNL1] = nonlinearTermsWithOldWVM(wvm,Ubar1,Vbar1,Wbar1,Nbar1);

[Ubar2,Vbar2,Wbar2,Nbar2] = physicalVarsFromWaveVortexVars(wvm,wvm.Ap,wvm.Am,wvm.A0);
[uNL2,vNL2,nNL2] = nonlinearTermsWithOldWVM(wvm,Ubar2,Vbar2,Wbar2,Nbar2);

% u_error = error3(U,U1);
% v_error = error3(V,V1);
% w_error = error3(W,W1);
% ux_error = error3(Ux,Ux1);
% vx_error = error3(Vx,Vx1);
% uy_error = error3(Uy,Uy1);
% vy_error = error3(Vy,Vy1);
% uz_error = error3(Uz,Uz1);
% vz_error = error3(Vz,Vz1);
% fprintf('\tu error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(u_error)))))));
% fprintf('\tv error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(v_error)))))));
% fprintf('\tw error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(w_error)))))));
% fprintf('\tux error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(ux_error)))))));
% fprintf('\tvx error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(vx_error)))))));
% fprintf('\tuy error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(uy_error)))))));
% fprintf('\tvy error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(vy_error)))))));
% fprintf('\tuz error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(uz_error)))))));
% fprintf('\tvz error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(vz_error)))))));

uNL_error = error3(uNL1,uNL2);
vNL_error = error3(vNL1,vNL2);
nNL_error = error3(nNL1,nNL2);

fprintf('\tuNL error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(uNL_error)))))));
fprintf('\tvNL error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(vNL_error)))))));
fprintf('\tnNL error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(nNL_error)))))));

[Fp1,Fm1,F01] = nonlinearFluxNewWVT(wvt,uNL1,vNL1,nNL1,wvm.EMAp);
[Fp2,Fm2,F02] = nonlinearFluxOldWVM(wvm,uNL2,vNL2,nNL2,wvm.EMAp);

Fp2_error = error3(Fp1,Fp2);
Fm2_error = error3(Fm1,Fm2);
F02_error = error3(F01,F02);

fprintf('\tFp error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(Fp2_error)))))));
fprintf('\tFm error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(Fm2_error)))))));
fprintf('\tF0 error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(F02_error)))))));

function [Ubar,Vbar,Wbar,Nbar] = physicalVarsFromWaveVortexVars(wvt,Ap,Am,A0)
Ubar = wvt.UAp.*Ap + wvt.UAm.*Am + wvt.UA0.*A0;
Vbar = wvt.VAp.*Ap + wvt.VAm.*Am + wvt.VA0.*A0;
Wbar = wvt.WAp.*Ap + wvt.WAm.*Am;
Nbar = wvt.NAp.*Ap + wvt.NAm.*Am + wvt.NA0.*A0;
end

function [uNL,vNL,nNL] = nonlinearTermsWithNewWVT(wvt,Ubar,Vbar,Wbar,Nbar)
[U,Ux,Uy,Uz] = wvt.transformToSpatialDomainWithFAllDerivatives(Ubar);
[V,Vx,Vy,Vz] = wvt.transformToSpatialDomainWithFAllDerivatives(Vbar);
W = wvt.transformToSpatialDomainWithG(Wbar);
[~,ETAx,ETAy,ETAz] = wvt.transformToSpatialDomainWithGAllDerivatives(Nbar);
uNL = -U.*Ux - V.*Uy - W.*Uz;
vNL = -U.*Vx - V.*Vy - W.*Vz;
nNL = -U.*ETAx - V.*ETAy - W.*ETAz;
end

function [uNL,vNL,nNL] = nonlinearTermsWithOldWVM(wvm,Ubar,Vbar,Wbar,Nbar)
[U,Ux,Uy,Uz] = wvm.TransformToSpatialDomainWithFAllDerivatives(Ubar);
[V,Vx,Vy,Vz] = wvm.TransformToSpatialDomainWithFAllDerivatives(Vbar);
W = wvm.TransformToSpatialDomainWithG(Wbar);
[~,ETAx,ETAy,ETAz] = wvm.TransformToSpatialDomainWithGAllDerivatives(Nbar);
uNL = -U.*Ux - V.*Uy - W.*Uz;
vNL = -U.*Vx - V.*Vy - W.*Vz;
nNL = -U.*ETAx - V.*ETAy - W.*ETAz;
end

function [Fp,Fm,F0] = nonlinearFluxOldWVM(wvm,uNL,vNL,nNL,AAMask)
uNLbar = wvm.TransformFromSpatialDomainWithF(uNL);
vNLbar = wvm.TransformFromSpatialDomainWithF(vNL);
nNLbar = wvm.TransformFromSpatialDomainWithG(nNL);

Fp = AAMask .* (wvm.ApU.*uNLbar + wvm.ApV.*vNLbar + wvm.ApN.*nNLbar);
Fm = AAMask .* (wvm.AmU.*uNLbar + wvm.AmV.*vNLbar + wvm.AmN.*nNLbar);
F0 = AAMask .* (wvm.A0U.*uNLbar + wvm.A0V.*vNLbar + wvm.A0N.*nNLbar);
end

function [Fp,Fm,F0] = nonlinearFluxNewWVT(wvt,uNL,vNL,nNL,AAMask)
uNLbar = wvt.transformFromSpatialDomainWithF(uNL);
vNLbar = wvt.transformFromSpatialDomainWithF(vNL);
nNLbar = wvt.transformFromSpatialDomainWithG(nNL);

Fp = AAMask .* (wvt.ApU.*uNLbar + wvt.ApV.*vNLbar + wvt.ApN.*nNLbar);
Fm = AAMask .* (wvt.AmU.*uNLbar + wvt.AmV.*vNLbar + wvt.AmN.*nNLbar);
F0 = AAMask .* (wvt.A0U.*uNLbar + wvt.A0V.*vNLbar + wvt.A0N.*nNLbar);
end
