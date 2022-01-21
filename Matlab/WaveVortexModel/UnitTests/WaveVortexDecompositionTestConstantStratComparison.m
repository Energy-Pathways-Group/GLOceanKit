%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% WaveVortexDecompositionTest
%
% This script tests the API decompose an existing (u,v,w,rho_prime) into
% wave-vortex components
%
% Jeffrey J. Early
% jeffrey@jeffreyearly.com
%
% April 12th, 2018      Version 1.0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 16;
aspectRatio = 1;

Lx = 100e3;
Ly = aspectRatio*100e3;
Lz = 5000;

Nx = N;
Ny = aspectRatio*N;
Nz = N+1; % 2^n + 1 grid points, to match the Winters model, but 2^n ok too.

latitude = 31;
N0 = 5.2e-3; % Choose your stratification 7.6001e-04

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wvm = WaveVortexModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0 );
oldwavemodel = InternalWaveModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Generate randomized amplitudes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ApIO,AmIO,ApIGW,AmIGW,A0G,A0G0,A0rhobar] = wvm.GenerateRandomFlowState();

Ap = ApIO + ApIGW;
Am = AmIO + AmIGW;
A0 = A0G + A0G0 + A0rhobar;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Forward/back transformation tests
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n********** Transform tests **********\n');

error = @(u,u_unit) max( [max(max(max(abs(u-u_unit)/max( [max(max(max( abs(u) ))), 1e-15] )))), 1e-15]);
error2 = @(u,u_unit) abs((u-u_unit))./(max(max(max(abs(u_unit)))));

t = 651;
[u,v,w,eta] = wvm.TransformWaveVortexToUVWEta(Ap,Am,A0,t);
[App,Amm,A00] = wvm.TransformUVEtaToWaveVortex(u,v,eta,t);

Ap_error = error2(Ap,App);
Am_error = error2(Am,Amm);
A0_error = error2(A0,A00);

fprintf('\tAp error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(Ap_error)))))));
fprintf('\tAm error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(Am_error)))))));
fprintf('\tA0 error: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(A0_error)))))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Forward/back transformation tests
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n********** Transform tests F/G **********\n');

% Having subsumed the coefficients for these transformations into the
% coefficients, these are no longer direct inverses. They should differ by
% a factor of 2*(Nz-1).
% First check the G transform
w_bar = wvm.TransformFromSpatialDomainWithG( w );
w_back = wvm.TransformToSpatialDomainWithG(w_bar);
w_error = error2(w,w_back);
fprintf('\tG-transform: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(w_error)))))));

% First check the F transform
u_bar = wvm.TransformFromSpatialDomainWithF( u );
u_back = wvm.TransformToSpatialDomainWithF(u_bar);
u_error = error2(u,u_back);
fprintf('\tF-transform: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(u_error)))))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Derivative test
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n********** Derivative tests **********\n');

[App,Amm,A00] = wvm.TransformUVEtaToWaveVortex(u,v,eta);
Ubar = wvm.UAp.*App + wvm.UAm.*Amm + wvm.UA0.*A00;

u_unit = wvm.TransformWaveVortexToUVWEta(App,Amm,A00);
ux_unit = DiffFourier(wvm.x,u_unit,1,1);
uy_unit = DiffFourier(wvm.y,u_unit,1,2);
uz_unit = DiffCosine(wvm.z,u_unit,1,3);

[u,ux,uy,uz] = wvm.TransformToSpatialDomainWithFAllDerivatives( Ubar );

u_error = error2(u,u_unit);
fprintf('\tNo-derivative: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(u_error)))))));
u_error = error2(ux,ux_unit);
fprintf('\tx-derivative: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(u_error)))))));
u_error = error2(uy,uy_unit);
fprintf('\ty-derivative: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(u_error)))))));
u_error = error2(uz,uz_unit);
fprintf('\tz-derivative: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(u_error)))))));

Nbar = wvm.NAp.*App + wvm.NAm.*Amm + wvm.NA0.*A00;

[~,~,~,eta_unit] = wvm.TransformWaveVortexToUVWEta(App,Amm,A00);
etax_unit = DiffFourier(wvm.x,eta_unit,1,1);
etay_unit = DiffFourier(wvm.y,eta_unit,1,2);
etaz_unit = DiffSine(wvm.z,eta_unit,1,3);

[eta,etax,etay,etaz] = wvm.TransformToSpatialDomainWithGAllDerivatives( Nbar );

u_error = error2(eta,eta_unit);
fprintf('\tNo-derivative: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(u_error)))))));
u_error = error2(etax,etax_unit);
fprintf('\tx-derivative: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(u_error)))))));
u_error = error2(etay,etay_unit);
fprintf('\ty-derivative: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(u_error)))))));
u_error = error2(etaz,etaz_unit);
fprintf('\tz-derivative: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(u_error)))))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Decomposition test
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n********** Decomposition tests **********\n');

error2 = @(u,u_unit) max(max(max( abs((u(abs(u_unit)>1e-15)-u_unit(abs(u_unit)>1e-15)))./abs(u_unit(abs(u_unit)>1e-15)) )));

% [ApIO,AmIO,ApIGW,AmIGW,A0G,A0G0,A0rhobar]
for i=1:11
    if i <= 5 % First we walk through the four types of solutions in isolation
        mask = zeros(5,1);
        mask(i)=1;
    elseif i <= 10  % now walk through triplets
        mask = ones(5,1);
        mask(i-5) = 0;
    else
        mask = ones(5,1);
    end
    
    wvm.Ap = mask(1)*ApIO + mask(2)*ApIGW;
    wvm.Am = mask(1)*AmIO + mask(2)*AmIGW;
    wvm.A0 = mask(3)*A0G + mask(4)*A0G0 + mask(5)*A0rhobar;
        
    fprintf('\nmask %d\n',i);

    fprintf('total integrated energy: %f m^3/s\n', wvm.totalEnergy);
    fprintf('total spectral energy: %f m^3/s\n', wvm.totalSpectralEnergy);
    
    [u,v,w,eta] = wvm.TransformWaveVortexToUVWEta(wvm.Ap,wvm.Am,wvm.A0,t);
    oldwavemodel.InitializeWithHorizontalVelocityAndIsopycnalDisplacementFields(t,u,v,eta);
    [u_wm,v_wm,w_wm,eta_wm] = oldwavemodel.VariableFieldsAtTime(t,'u','v','w','zeta');
    u_error = error(u,u_wm);
    fprintf('\n\tu: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(u_error)))))));
    u_error = error(v,v_wm);
    fprintf('\tv: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(u_error)))))));
    u_error = error(w,w_wm);
    fprintf('\tw: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(u_error)))))));
    u_error = error(eta,eta_wm);
    fprintf('\teta: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(u_error)))))));


    fprintf('\n\n');
end