%% 2024-10-17
% Informal unit test to confirm that we are able to correctly create the
% vertical mode matrices.

Lxyz = [1000, 500, 500];
Nxyz = [16 8 9];
latitude = 33;
wvt_c = WVTransformConstantStratification(Lxyz, Nxyz, latitude=latitude, isHydrostatic=0);
wvt_b = WVTransformBoussinesq(Lxyz, Nxyz, N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));

%%
Finv_c = wvt_c.FinvMatrix;
Finv_b = wvt_b.FinvMatrix;

Ginv_c = wvt_c.GinvMatrix;
Ginv_b = wvt_b.GinvMatrix;

max(abs(Finv_b(:)-Finv_c(:)))

%% Go compute the vertical modes with the correct normalization
verticalModes = InternalModesWKBSpectral(N2=wvt_b.N2Function,zIn=[-wvt_b.Lz 0],zOut=wvt_b.z,latitude=wvt_b.latitude,nModes=wvt_b.Nj-1,nEVP=max(256,floor(2.1*wvt_b.Nz)));
verticalModes.normalization = Normalization.geostrophic;
verticalModes.upperBoundary = UpperBoundary.rigidLid;
[Finv,Ginv,h] = verticalModes.ModesAtFrequency(0);

h_a = cat(2,1,h);
Finv_a = cat(2,ones(wvt_b.Nz,1),Finv);
Ginv_a = cat(2,zeros(wvt_b.Nz,1),Ginv);

max(abs(Finv_b(:)-Finv_c(:)))
max(abs(Ginv_b(:)-Ginv_a(:)))
max(abs(Ginv_c(:)-Ginv_a(:)))
max(abs(wvt_b.h_0(:)-h_a(:)))
max(abs(wvt_c.h_0(:)-h_a(:)))

%%

iK = 15;
k = wvt_c.Kh(1,iK);
[kMode,lMode,j] = wvt_b.modeNumberFromIndex(iK);
verticalModes.normalization = Normalization.kConstant;
verticalModes.upperBoundary = UpperBoundary.rigidLid;
[Finv,Ginv,h] = verticalModes.ModesAtWavenumber(k);

h_a = cat(2,1,h);
FwInv_a = cat(2,ones(wvt_b.Nz,1),Finv);
GwInv_a = cat(2,zeros(wvt_b.Nz,1),Ginv);

%%
[kMode,lMode] = wvt_c.horizontalModes.modeNumberFromIndex(iK);
FwInv_b = wvt_b.FwInvMatrix(kMode,lMode);
FwInv_c = wvt_c.FwInvMatrix(kMode,lMode);
Fw_b = wvt_b.FwMatrix(kMode,lMode);
Fw_c = wvt_c.FwMatrix(kMode,lMode);
GwInv_b = wvt_b.GwInvMatrix(kMode,lMode);
GwInv_c = wvt_c.GwInvMatrix(kMode,lMode);
Gw_b = wvt_b.GwMatrix(kMode,lMode);
Gw_c = wvt_c.GwMatrix(kMode,lMode);
h_b = wvt_b.h_pm(:,iK);
h_c = wvt_c.h_pm(:,iK);

max(abs(FwInv_b(:)-FwInv_a(:)))
max(abs(FwInv_c(:)-FwInv_a(:)))
max(abs(GwInv_b(:)-GwInv_a(:)))
max(abs(GwInv_c(:)-GwInv_a(:)))
max(abs(h_b(:)-h_a(:)))
max(abs(h_b(:)-h_a(:)))

max(max(abs(Fw_b * FwInv_b - eye(wvt_b.Nj))))
max(max(abs(Fw_c * FwInv_c - eye(wvt_c.Nj))))
A = Gw_c * GwInv_c;
max(max(abs(A(2:end,2:end) - eye(wvt_c.Nj-1))))
A = Gw_b * GwInv_b;
max(max(abs(A(2:end,2:end) - eye(wvt_c.Nj-1))))

%%

wvt_c.transformWithG_wg(ones(wvt_c.spectralMatrixSize))
wvt_b.transformWithG_wg(ones(wvt_c.spectralMatrixSize))