L = 100;
Lxyz = [2*pi*L, 2*pi*L, pi*L];
Nxyz = [8 8 9];
latitude = 33;

wvt = WVTransformConstantStratification(Lxyz, Nxyz, latitude=latitude, isHydrostatic=0,shouldAntialias=0);

% insert the triad
wvt.initWithWaveModes(kMode=1,lMode=0,j=1,phi=0,u=0.05,sign=1); % wave A
wvt.addWaveModes(kMode=1,lMode=1,j=1,phi=0,u=0.05,sign=1); % wave B
wvt.addWaveModes(kMode=0,lMode=1,j=2,phi=0,u=0.05,sign=1); % wave C

indexA = wvt.indexFromModeNumber(1,0,1); % (2,3)
indexB = wvt.indexFromModeNumber(1,1,1); % (2,5)
indexC = wvt.indexFromModeNumber(0,1,2); % (3,2)

gamma = (wvt.N0^2 - wvt.f^2)/wvt.g;
Aw = sqrt(2/(wvt.Lz*gamma));

% non-dimensional periods of each wave
Ta = wvt.f/wvt.Omega(indexA);
Tb = wvt.f/wvt.Omega(indexB);
Tc = wvt.f/wvt.Omega(indexC);

% equivalent depth
ha = wvt.h_pm(indexA);
hb = wvt.h_pm(indexB);
hc = wvt.h_pm(indexC);

% vertical wave number
ma = wvt.J(indexA)*pi/wvt.Lz;
mb = wvt.J(indexB)*pi/wvt.Lz;
mc = wvt.J(indexC)*pi/wvt.Lz;

% The true amplitude has a factor of 2 (half-complex)
A = 2*wvt.Ap(indexA)*ha*Aw;
B = 2*wvt.Ap(indexB)*hb*Aw/sqrt(2);
C = 2*wvt.Ap(indexC)*hc*Aw;

% Define abs and relative error metrics
error = @(x,y) max(abs(x(:)-y(:)));
error = @(x,y) fprintf('relative error: %g\n',max(abs(x(:)-y(:))./abs(x(:))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Check that our analytical expressions for (u,v,w,eta) match
% These are pulled from the analytical expressions in my notes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u_bar = wvt.DCT * wvt.transformFromSpatialDomainWithFourier(wvt.u);
v_bar = wvt.DCT * wvt.transformFromSpatialDomainWithFourier(wvt.v);
w_bar = wvt.DST * wvt.transformFromSpatialDomainWithFourier(wvt.w);
n_bar = wvt.DST * wvt.transformFromSpatialDomainWithFourier(wvt.eta);

error( u_bar(indexA), (A/L/2)*(-1))
error( u_bar(indexB), (B/L/2)*(-1 + sqrt(-1)*Tb))
error( u_bar(indexC), (C/L/2)*(-sqrt(-1)*2*Tc))

error( v_bar(indexA), (A/L/2)*(-sqrt(-1)*Ta))
error( v_bar(indexB), (B/L/2)*(-1 - sqrt(-1)*Tb))
error( v_bar(indexC), (C/L/2)*(2))

error( w_bar(indexA), (A/L/2)*sqrt(-1))
error( w_bar(indexB), (B/L/2)*2*sqrt(-1))
error( w_bar(indexC), (C/L/2)*(-sqrt(-1)))

error( n_bar(indexA), (A/L/2/wvt.f)*Ta)
error( n_bar(indexB), (B/L/2/wvt.f)*2*Tb)
error( n_bar(indexC), (C/L/2/wvt.f)*(-Tc))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Check that our analytical expressions for (uNl,vNL,wNL,etaNL) match
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uNL = wvt.u .* wvt.diffX(wvt.u)   + wvt.v .* wvt.diffY(wvt.u)   + wvt.w .*  wvt.diffZF(wvt.u);
vNL = wvt.u .* wvt.diffX(wvt.v)   + wvt.v .* wvt.diffY(wvt.v)   + wvt.w .*  wvt.diffZF(wvt.v);
wNL = wvt.u .* wvt.diffX(wvt.w)   + wvt.v .* wvt.diffY(wvt.w)   + wvt.w .*  wvt.diffZG(wvt.w);
nNL = wvt.u .* wvt.diffX(wvt.eta) + wvt.v .* wvt.diffY(wvt.eta) + wvt.w .*  wvt.diffZG(wvt.eta);
uNL_bar = wvt.DCT * wvt.transformFromSpatialDomainWithFourier(uNL);
vNL_bar = wvt.DCT * wvt.transformFromSpatialDomainWithFourier(vNL);
wNL_bar = wvt.DST * wvt.transformFromSpatialDomainWithFourier(wNL);
nNL_bar = wvt.DST * wvt.transformFromSpatialDomainWithFourier(nNL);

% Renormalization: Must divide each amplitude by L AND because of the
% derivative also divide by L
uNL_A = (B*C/8/L/L/L)*(8*Tc-Tb-sqrt(-1)*(4*Tb*Tc+1));
uNL_B = (A*C/8/L/L/L)*(-6*Tc-sqrt(-1)*(1+2*Ta*Tc));
uNL_C = (A*B/8/L/L/L)*(Ta+Tb-sqrt(-1)*(1+Ta*Tb));

error(uNL_bar(indexA),uNL_A)
error(uNL_bar(indexB),uNL_B)
error(uNL_bar(indexC),uNL_C)

vNL_A = (B*C/8/L/L/L)*(2*Tc - Tb + sqrt(-1)*(2*Tb*Tc-7));
vNL_B = (A*C/8/L/L/L)*(3*Ta + sqrt(-1)*(-2*Ta*Tc-4));
vNL_C = (A*B/8/L/L/L)*(-2*Ta-2*Tb + sqrt(-1)*(2*Ta*Tb+2));
error(vNL_bar(indexA),vNL_A)
error(vNL_bar(indexB),vNL_B)
error(vNL_bar(indexC),vNL_C)

wNL_A = (B*C/8/L/L/L)*(5 + sqrt(-1)*(-Tb + 4*Tc));
wNL_B = (A*C/8/L/L/L)*(-1 + sqrt(-1)*(-Ta - 2*Tc));
wNL_C = (A*B/8/L/L/L)*(7 + sqrt(-1)*(-2*Ta-Tb));
error(wNL_bar(indexA),wNL_A)
error(wNL_bar(indexB),wNL_B)
error(wNL_bar(indexC),wNL_C)

nNL_A = (B*C/8/L/L/L/wvt.f)*(5*Tb*Tc + sqrt(-1)*(-2*Tb + 3*Tc));
nNL_B = (A*C/8/L/L/L/wvt.f)*(-3*Ta*Tc + sqrt(-1)*(-Ta + 2*Tc));
nNL_C = (A*B/8/L/L/L/wvt.f)*(-Ta*Tb + sqrt(-1)*(3*Ta-4*Tb));
error(nNL_bar(indexA),nNL_A)
error(nNL_bar(indexB),nNL_B)
error(nNL_bar(indexC),nNL_C)