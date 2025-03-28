Lxyz = [1000, 500, 500];
Nxyz = [32 16 17];
wvt = WVTransformStratifiedQG(Lxyz, Nxyz, N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));

% %%
% u = wvt.u;
% return
%%
[X,Y,Z] = wvt.xyzGrid;
Lx = wvt.Lx;
Ly = wvt.Ly;
phix = 2*pi*rand(1);
phiy = 2*pi*rand(1);

k_n = 1; l_n = 1;
kx = 2*pi*k_n/Lx;
ky = 2*pi*l_n/Ly;
f = cos(kx*X+phix) .* cos(ky*Y+phiy);

Df_analytical = -kx*sin(kx*X+phix).*cos(ky*Y+phiy);
Df_numerical = wvt.diffX(f);
max(abs(Df_analytical(:)-Df_numerical(:)))

Df_analytical = -ky*sin(ky*Y+phiy).*cos(kx*X+phix);
Df_numerical = wvt.diffY(f);
max(abs(Df_analytical(:)-Df_numerical(:)))