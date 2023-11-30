% This test suite tests the following methods:
% -diffX
% -diffY
% -diffZF
% -diffZG
% We do not compute ANY Nyquist modes, since they are not fully resolved.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Fourier derivatives
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Testing Fourier transform:\n')

Nx = 16;
Ny = 8;
Nz = 16+1;

Lx = 1.0;
Ly = 10.0;
Lz = 4.0;

latitude = 25;
N0 = 5.2e-3; % Choose your stratification 7.6001e-04
wvt = WVTransformConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], N0, latitude=latitude);

[X,Y,Z] = wvt.xyzGrid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transformWithF, 1st-dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi = pi/8;
f = @(k) cos(k*X+phi) .* cos(0*Y) .* cos(0*Z);
iK = 2*pi*(1:floor(Nx/2))'/Lx;

Df_numerical = @(u) wvt.transformToSpatialDomainWithF(A0=wvt.transformFromSpatialDomainWithF(u));
testname = sprintf('transformWithF, 1st dimension');
ReportErrors(f,f,Df_numerical,testname,iK(1:end-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transformWithF, 2nd-dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi = pi/8;
f = @(l) cos(0*X) .* cos(l*Y+phi) .* cos(0*Z);
iL = 2*pi*(1:floor(Ny/2))'/Ly;

Df_numerical = @(u) wvt.transformToSpatialDomainWithF(A0=wvt.transformFromSpatialDomainWithF(u));
testname = sprintf('transformWithF, 2nd dimension');
ReportErrors(f,f,Df_numerical,testname,iL(1:end-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transformWithF, 3rd-dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = @(m) cos(0*X) .* cos(0*Y) .* cos(m*Z);
dm = 1/((Nz-1)*(wvt.z(2)-wvt.z(1)));
iM = pi*dm*(0:(Nz-1))';

Df_numerical = @(u) wvt.transformToSpatialDomainWithF(A0=wvt.transformFromSpatialDomainWithF(u));
testname = sprintf('transformWithF, 3rd dimension');
ReportErrors(f,f,Df_numerical,testname,iM(1:end-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transformWithF-diffX, 1st-dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi = pi/8;
f = @(k) cos(k*X+phi) .* cos(0*Y) .* cos(0*Z);
Df_analytical = @(k) -k*sin(k*X+phi).* cos(0*Y) .* cos(0*Z);
iK = 2*pi*(1:floor(Nx/2))'/Lx;

Df_numerical = @(u) diffXF(wvt,u);
testname = sprintf('transformWithF-diffX, 1st dimension');
ReportErrors(f,Df_analytical,Df_numerical,testname,iK(1:end-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transformWithF-diffY, 2nd-dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi = pi/8;
f = @(l) cos(0*X) .* cos(l*Y+phi) .* cos(0*Z);
Df_analytical = @(l) -l*cos(0*X) .* sin(l*Y+phi) .* cos(0*Z);
iL = 2*pi*(1:floor(Ny/2))'/Ly;

Df_numerical = @(u) diffYF(wvt,u);
testname = sprintf('transformWithF-diffY, 2nd dimension');
ReportErrors(f,Df_analytical,Df_numerical,testname,iL(1:end-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transformWithF-diffZ, 3rd-dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = @(m) cos(0*X) .* cos(0*Y) .* cos(m*Z);
Df_analytical = @(m) -m*cos(0*X) .* cos(0*Y) .* sin(m*Z);
dm = 1/((Nz-1)*(wvt.z(2)-wvt.z(1)));
iM = pi*dm*(0:(Nz-1))';

Df_numerical = @(u) diffZF(wvt,u);
testname = sprintf('transformWithF-diffZ, 3rd dimension');
ReportErrors(f,Df_analytical,Df_numerical,testname,iM(1:end-1));

function ux = diffXF(wvt,u)
[u,ux,uy,uz] = wvt.transformToSpatialDomainWithFAllDerivatives_FFT(wvt.transformFromSpatialDomainWithF_FFT(u));
end

function uy = diffYF(wvt,u)
[u,ux,uy,uz] = wvt.transformToSpatialDomainWithFAllDerivatives_FFT(wvt.transformFromSpatialDomainWithF_FFT(u));
end

function uz = diffZF(wvt,u)
[u,ux,uy,uz] = wvt.transformToSpatialDomainWithFAllDerivatives_FFT(wvt.transformFromSpatialDomainWithF_FFT(u));
end





