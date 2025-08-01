N0 = 5.2e-3;
wvtc = WVTransformConstantStratification([4e3, 4e3, 2e3], [16 16 9], N0=N0,shouldAntialias=false);
wvt = WVTransformBoussinesq([4e3, 4e3, 2e3], [16 16 9], N2=@(z)N0*N0*ones(size(z)),shouldAntialias=false);

wvtc.initWithWaveModes(kMode=1,lMode=0,j=1,phi=0,u=0.05,sign=1); % wave A
wvtc.addWaveModes(kMode=1,lMode=1,j=1,phi=0,u=0.05,sign=1); % wave B
wvtc.addWaveModes(kMode=0,lMode=1,j=2,phi=0,u=0.05,sign=1); % wave C

wvt.initWithWaveModes(kMode=1,lMode=0,j=1,phi=0,u=0.05,sign=1); % wave A
wvt.addWaveModes(kMode=1,lMode=1,j=1,phi=0,u=0.05,sign=1); % wave B
wvt.addWaveModes(kMode=0,lMode=1,j=2,phi=0,u=0.05,sign=1); % wave C

%%
[Fp,Fm,F0] = wvtc.nonlinearFlux();

%%
Ac = wvtc.transformWithG_wg(ones(wvtc.spectralMatrixSize));
A = wvt.transformWithG_wg(ones(wvtc.spectralMatrixSize));