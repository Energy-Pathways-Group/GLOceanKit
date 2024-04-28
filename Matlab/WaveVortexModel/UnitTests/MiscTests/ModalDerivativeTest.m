Lxy = 750e3;
N = 16;
Nz = 60;

Lz = 4000;
N0 = 3*2*pi/3600; % buoyancy frequency at the surface, radians/seconds
L_gm = 1300; % thermocline exponential scale, meters
N2 = @(z) N0*N0*exp(2*z/L_gm);
wvt = WVTransformHydrostatic([Lxy, Lxy, Lz], [N, N, Nz], N2=N2,latitude=31);

% coeff = (wvt.Q ./ wvt.P);

%%
Finv = wvt.FinvMatrix;
f = Finv(:,3)/max(Finv(:,3));
figure, plot(f);

a = shiftdim(f,-2) .* ones(size(wvt.X));
da = wvt.diffZF(a);

figure, plot(squeeze(da(1,1,:)),wvt.z)
hold on, plot(diff(squeeze(a(1,1,:)))./diff(wvt.z),wvt.z(2:end))