% wvt = WVTransformBoussinesq([15e3, 15e3, 5000], [16 16 5], N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));
% wvt = WVTransformConstantStratification([15e3, 15e3, 5000], [8 8 5]);
% wvt = WVTransformHydrostatic([15e3, 15e3, 5000], [8 8 5], N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));

wvt = WVTransformConstantStratification([15e3, 15e3, 1300], [32 32 17]);

% try
%     wvt = WVTransformConstantStratification([15e3, 15e3, 5000], [8 8 5],latitude=90);
% catch ME
%     disp(ME)
% end

%%
% flowComponent = WVMeanDensityAnomalyComponent(wvt);
flowComponent = WVInternalGravityWaveComponent(wvt);
% flowComponent = WVInertialOscillationComponent(wvt);
% flowComponent = WVGeostrophicComponent(wvt);
[Ap,Am,A0] = flowComponent.randomAmplitudes;
varianceMatrix = abs(Ap).^2 + abs(Am).^2 + abs(A0).^2;
radialVarianceMatrix = wvt.transformToRadialWavenumber(varianceMatrix);
sum(radialVarianceMatrix(:))
sum(varianceMatrix(:))

%%
solutionIndex = 1; %% 2, 6, 14, 24, 40, 46, 58, 68, lMode=0!!!
soln = flowComponent.solutionForModeAtIndex(solutionIndex,amplitude='random');

wvt.t = 6000;
args = {wvt.X,wvt.Y,wvt.Z,wvt.t};
wvt.initWithUVEta(soln.u(args{:}), soln.v(args{:}),soln.eta(args{:}));

wvt.t = 86400;
args = {wvt.X,wvt.Y,wvt.Z,wvt.t};
% u = wvt.u;
% u2 = soln.u(args{:});
% v = wvt.v;
% v2 = soln.v(args{:});
% p = wvt.p;
% p2 = soln.p(args{:});
% w = wvt.w;
% w2 = soln.w(args{:});
% eta = wvt.eta;
% eta2 = soln.eta(args{:});
qgpv = wvt.qgpv;
qgpv2 = soln.qgpv(args{:});

a = soln.depthIntegratedTotalEnergy(isHydrostatic=wvt.isHydrostatic);
b = soln.energyFactor*(abs(soln.coefficientMatrixAmplitude).^2 + abs(soln.conjugateCoefficientMatrixAmplitude).^2);