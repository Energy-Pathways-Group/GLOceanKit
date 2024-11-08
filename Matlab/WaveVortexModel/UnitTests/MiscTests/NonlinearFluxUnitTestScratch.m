L = 100;
Lxyz = [2*pi*L, 2*pi*L, pi*L];
Nxyz = [8 8 9];
latitude = 33;

wvt = WVTransformConstantStratification(Lxyz, Nxyz, latitude=latitude, isHydrostatic=0,shouldAntialias=0);

% insert the triad
wvt.initWithWaveModes(kMode=1,lMode=0,j=1,phi=0,u=0.05,sign=1);
wvt.addWaveModes(kMode=1,lMode=1,j=1,phi=0,u=0.05,sign=1);
wvt.addWaveModes(kMode=0,lMode=1,j=2,phi=0,u=0.05,sign=1);

index101 = wvt.indexFromModeNumber(1,0,1); % (2,3)
index111 = wvt.indexFromModeNumber(1,1,1); % (2,5)
index012 = wvt.indexFromModeNumber(0,1,2); % (3,2)

gamma = (wvt.N0^2 - wvt.f^2)/wvt.g;
Aw = sqrt(2/(wvt.Lz*gamma));

Ta = wvt.f/wvt.Omega(index101);
Tb = wvt.f/wvt.Omega(index111);
Tc = wvt.f/wvt.Omega(index012);

ha = wvt.h_pm(index101);
hb = wvt.h_pm(index111);
hc = wvt.h_pm(index012);

ma = wvt.J(index101)*pi/wvt.Lz;
mb = wvt.J(index111)*pi/wvt.Lz;
mc = wvt.J(index012)*pi/wvt.Lz;

% The true amplitude has a factor of 2 (half complex)
A = 2*wvt.Ap(index101)*ha*Aw;
B = 2*wvt.Ap(index111)*hb*Aw/sqrt(2);
C = 2*wvt.Ap(index012)*hc*Aw;

uNL = wvt.u .* wvt.diffX(wvt.u)   + wvt.v .* wvt.diffY(wvt.u)   + wvt.w .*  wvt.diffZF(wvt.u);
uNL_bar = wvt.DCT * wvt.transformFromSpatialDomainWithFourier(uNL);

AuNL = uNL_bar(index101);
BuNL = uNL_bar(index111);
CuNL = uNL_bar(index012);

% Renormalization: Must divide each amplitude by L AND because of the
% derivative also divide by L
AuNL_expected = -B*C*(Tb-8*Tc+sqrt(-1)*(4*Tb*Tc+1))/8/L/L/L;
BuNL_expected = -A*C*(6*Tc+sqrt(-1)*(1+2*Ta*Tc))/8/L/L/L;
CuNL_expected = A*B*(Ta+Tb-sqrt(-1)*(1+Ta*Tb))/8/L/L/L;

u_bar = wvt.DCT * wvt.transformFromSpatialDomainWithFourier(wvt.u);
v_bar = wvt.DCT * wvt.transformFromSpatialDomainWithFourier(wvt.v);
w_bar = wvt.DST * wvt.transformFromSpatialDomainWithFourier(wvt.w);

error = @(x,y) max(abs(x(:)-y(:)));
error = @(x,y) max(abs(x(:)-y(:))./abs(x(:)));

error( u_bar(index101), -A/L/2)
error( u_bar(index111), -B/L/2 + sqrt(-1)*B*Tb/L/2)
error( u_bar(index012), -sqrt(-1)*2*C*Tc/L/2)

error( v_bar(index101), -sqrt(-1)*A*Ta/L/2)
error( v_bar(index111), B*(-1 - sqrt(-1)*Tb)/L/2)
error( v_bar(index012), 2*C/L/2)

error(w_bar(index101), A*sqrt(-1)/L/2)
error(w_bar(index111), B*2*sqrt(-1)/L/2)
error(w_bar(index012), -C*sqrt(-1)/L/2)

error(AuNL,AuNL_expected)
error(BuNL,BuNL_expected)
error(CuNL,CuNL_expected) % This on is not okay!