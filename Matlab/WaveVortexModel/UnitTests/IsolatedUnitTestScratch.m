wvt = WVTransformHydrostatic([15e3, 15e3, 5000], [8 8 5], N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));
solutionGroup = WVGeostrophicSolutionGroup(wvt);
args = {wvt.X,wvt.Y,wvt.Z,wvt.t};

%%
solutionIndex = 1;
soln = solutionGroup.uniqueSolutionAtIndex(solutionIndex,amplitude='random');

wvt.initWithUVEta(soln.u(args{:}), soln.v(args{:}),soln.eta(args{:}));
u = wvt.u;
u2 = soln.u(args{:});
v = wvt.v;
v2 = soln.v(args{:});