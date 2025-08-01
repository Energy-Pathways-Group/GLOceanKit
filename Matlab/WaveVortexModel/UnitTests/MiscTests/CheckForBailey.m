file = '/Volumes/Samsung_T5/nsf_iwv/WintersNonlinear/EarlyV2_GM_NL_forced_damped_01xGM';
WM = WintersModel(file);

[x,y,z,rho_bar] = WM.VariableFieldsFrom3DOutputFileAtIndex(1,'x','y','z','s1_bar');
Nx = length(x);
Ny = length(y);
Nz = length(z);
dx = (max(x)-min(x))/(Nx-1);
Lx = dx*Nx;
dy = (max(y)-min(y))/(Ny-1);
Ly = dy*Ny;
Lz = max(z)-min(z);
f = WM.VariableFieldsFrom3DOutputFileAtIndex(1,'f');
N0 = InternalModesConstantStratification.BuoyancyFrequencyFromConstantStratification(rho_bar,z);
latitude = double(asind(f/(2 * 7.2921E-5)));

wavemodel = InternalWaveModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0, min(rho_bar));
wvm = WaveVortexModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0, min(rho_bar));

[t,u,v,w,rho_prime] = WM.VariableFieldsFrom3DOutputFileAtIndex(fileIncrements(iTime),'t','u','v','w','rho_prime');

wavemodel.InitializeWithHorizontalVelocityAndDensityPerturbationFields(t,u,v,rho_prime);