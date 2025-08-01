 % wvt = WVTransformBoussinesq([15e3, 15e3, 5000], [16 16 5], N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));
% wvt = WVTransformConstantStratification([15e3, 15e3, 5000], [8 8 5]);
% wvt = WVTransformHydrostatic([15e3, 15e3, 5000], [32 32 10], N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));

wvt = WVTransformConstantStratification([15e3, 15e3, 1300], [16 16 9]);

% try
%     wvt = WVTransformConstantStratification([15e3, 15e3, 5000], [8 8 5],latitude=90);
% catch ME
%     disp(ME)
% end
%%
solutionIndex = 1;

args = {wvt.X,wvt.Y,wvt.Z,wvt.t};
for solutionIndex = 1:1 %wvt.flowComponent("wave").nModes
    soln = wvt.flowComponent("wave").solutionForModeAtIndex(solutionIndex,amplitude='random');
    wvt.initWithUVEta(soln.u(args{:}), soln.v(args{:}),soln.eta(args{:}));
    [Fp,Fm,F0] = wvt.nonlinearFlux();
    [Ep,Em,E0] = wvt.energyFluxFromNonlinearFlux(Fp,Fm,F0,deltaT=10);
    fprintf('(Ep,Em,E0) = (%g, %g, %g), total %g\n', sum(Ep(:)), sum(Em(:)), sum(E0(:)), sum(Ep(:))+sum(Em(:))+sum(E0(:)));
end

%%
U_io = 0.2;
Ld = wvt.Lz/2;
theta = 0;
u_NIO = @(z) U_io*cos(theta)*exp((z/Ld));
v_NIO = @(z) U_io*sin(theta)*exp((z/Ld));

wvt.initWithInertialMotions(u_NIO,v_NIO);
u1 = wvt.u_io;
v1 = wvt.v_io;
Ap1 = wvt.Ap; Am1 = wvt.Am; A01 = wvt.A0;

%% Populate the flow field with junk...
wvt.initWithRandomFlow();
u2 = wvt.u_io;
v2 = wvt.v_io;
Ap2 = wvt.Ap; Am2 = wvt.Am; A02 = wvt.A0;

%%
wvt.addInertialMotions(u_NIO,v_NIO);
u3 = wvt.u_io;
v3 = wvt.v_io;
Ap3 = wvt.Ap; Am3 = wvt.Am; A03 = wvt.A0;


%%
self.verifyThat(u1 + u2,IsSameSolutionAs(u3),'u_tot');
self.verifyThat(v1 + v2,IsSameSolutionAs(v3),'v_tot');

%%
wvt.removeAll();
kMode = -3; lMode = -5; j=3; phi = pi*0.3; u=0.1; sign = +1;
wvt.t = 22654;
wvt.setWaveModes(kMode=kMode,lMode=lMode,j=j,phi=phi,u=u,sign=sign);
soln = wvt.waveComponent.internalGravityWaveSolution(kMode,lMode,j,u,phi,sign,amplitudeIsMaxU=1,t=wvt.t);
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