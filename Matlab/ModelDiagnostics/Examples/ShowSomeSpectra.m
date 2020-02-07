file = '/Volumes/Samsung_T5/nsf_iwv/WintersNonlinear/EarlyV2_GM_NL_forced_damped_restart';

WM=WintersModel(file);
[t0,u,v,w,rho_prime] = WM.VariableFieldsFrom3DOutputFileAtIndex(1,'t','u','v','w','rho_prime');

% extract the variables, we will automate this later.
wavemodel = WM.wavemodel;
dims = [wavemodel.Lx,wavemodel.Ly,wavemodel.Lz];
n = [wavemodel.Nx,wavemodel.Ny,wavemodel.Nz];
latitude = wavemodel.latitude;
rhoBar = wavemodel.RhoBarAtDepth(wavemodel.z);
N2 = wavemodel.N2AtDepth(wavemodel.z);

diag = ModelDiagnostics(dims, n, latitude, rhoBar, N2);
diag.InitializeWithHorizontalVelocityAndDensityPerturbationFields(u,v,w,rho_prime);

S = diag.MakeIsotropicWavenumberSpectrumFromXY(u+sqrt(-1)*v);

SsumVsDepth = sum(S,1)*diag.dK;
[Euv,z] = diag.HorizontalVelocityVariance();

figure, plot(1e4*SsumVsDepth,z), hold on, plot(1e4*Euv,z)
xlabel('cm^2/s^2')

