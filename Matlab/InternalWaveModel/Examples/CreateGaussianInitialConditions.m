Lx = 15e3;
Ly = 15e3;
Lz = 5000;

Nx = 6;
Ny = 8;
Nz = 4;

latitude = 31;
N0 = 5.2e-3/2; % Choose your stratification 7.6001e-04
U = 0.01; % m/s
phi = 0*0.232; % random, just to see if its doing the right thing
t = 0*86400;

wavemodel = InternalWaveModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);

k_loop = 2;
l_loop = 2;
j0 = 3;
sign = 1;
wavemodel.InitializeWithPlaneWave(k_loop,l_loop,j0,U,sign);

zeta = wavemodel.IsopycnalDisplacementFieldAtTime(t);
zeta_bar = wavemodel.TransformFromSpatialDomainWithG(zeta);
zeta_bar(abs(zeta_bar)/(max(max(max(abs(zeta_bar))))) < 1e-7) = 0;

Kh = wavemodel.Kh;
Kh(Kh<1e-14) = 1;
negate_zeta = abs(wavemodel.Omega)./(Kh .* sqrt(wavemodel.h));
A_plus = -zeta_bar .* negate_zeta ./ wavemodel.G; % extra factor from constant stratification case.
