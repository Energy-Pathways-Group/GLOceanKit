L = 100;
Lxyz = [2*pi*L, 2*pi*L, pi*L];
Nxyz = [8 8 9];
latitude = 33;

wvt = WVTransformConstantStratification(Lxyz, Nxyz, latitude=latitude, isHydrostatic=0,shouldAntialias=0);

k=1; l=0; j=1;
wvt.addWaveModes(kMode=k,lMode=l,j=j,phi=0,u=0.05,sign=1);
index = wvt.indexFromModeNumber(k,l,j);
wvt.Ap(index)


T = sqrt( (k^2 + l^2 + j^2)/( ((wvt.N0/wvt.f)^2)*(k^2 + l^2) + j^2 ));

u = @(x,y,z) cos(k*x/(2*pi*L))