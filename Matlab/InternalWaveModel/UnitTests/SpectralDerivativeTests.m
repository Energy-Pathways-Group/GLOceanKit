Nx = 16;
Ny = 8;
Nz = 9;

Lx = 10;
Ly = 5;
Lz = 5;

wavemodel = InternalWaveModelConstantStratification([Lx Ly Lz], [Nx Ny Nz], 33, 5.2e-3);

X = wavemodel.X;

f = 1*ones(size(X)) + 3*sin(2*(2*pi)*X);
[f_bar, K, L, M] = wavemodel.TransformFromSpatialDomainWithFFull(f);
fx_bar = K .* f_bar;
f_x = wavemodel.TransformToSpatialDomainWithFFull(fx_bar);