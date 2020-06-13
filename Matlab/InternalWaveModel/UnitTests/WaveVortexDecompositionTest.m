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

rho0 = 1025; g = 9.81;
rho = @(z) -(N0*N0*rho0/g)*z + rho0;
z = linspace(-Lz,0,Nz);
shouldUseArbitraryStratificationModel = 1;

rho = InternalModes.StratificationProfileWithName('exponential');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if shouldUseArbitraryStratificationModel == 0
    wavemodel = InternalWaveModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);
else
    wavemodel = InternalWaveModelArbitraryStratification([Lx, Ly, Lz], [Nx, Ny, Nz], rho, z, latitude, 'method','spectral','cacheFile','evp-cache.mat'); %, 'method','spectral'
end

% wavemodel.InitializeWithGMSpectrum(1.0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Generate randomized amplitudes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We have to exclude the nyquist, because it's no resolvable.
shouldExcludeNyquist = 1;

A0 = randn(1,1)+sqrt(-1)*randn(1,1);

% Generate waves---with one additional condition for the inertial waves.
Ap = InternalWaveModel.GenerateHermitianRandomMatrix( size(wavemodel.K), shouldExcludeNyquist );
Am = InternalWaveModel.GenerateHermitianRandomMatrix( size(wavemodel.K), shouldExcludeNyquist );
Am(1,1,:) = conj(Ap(1,1,:)); % Inertial motions go only one direction!

% Generate barotropic geostrophic currents--amplitude prefactor of 1e-3 is
% set in order to keep the energetics similar magnitude to the waves.
B0 = 2e-3*InternalWaveModel.GenerateHermitianRandomMatrix( size(wavemodel.K(:,:,1)), shouldExcludeNyquist);
B0(1,1) = 0;

% Generate internal geostrophic currents--amplitude prefactor of 6e-2 is
% set in order to keep the energetics similar magnitude to the waves.
B = 6e-2*InternalWaveModel.GenerateHermitianRandomMatrix( size(wavemodel.K), shouldExcludeNyquist );
B(1,1,:) = 0;

% Now initial th models with these.
wavemodel.A0 = A0;
wavemodel.GenerateWavePhases(Ap,Am);
wavemodel.GenerateGeostrophicCurrents(B0,B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Forward/back transformation tests
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n********** Transform tests **********\n');

error = @(u,u_unit) max( [max(max(max(abs(u-u_unit)/max( [max(max(max( abs(u) ))), 1e-15] )))), 1e-15]);
error2 = @(u,u_unit) abs((u-u_unit))./(max(max(max(abs(u_unit)))));

% First check the baroclinic G transform
w = wavemodel.TransformToSpatialDomainWithG( wavemodel.w_plus );
w_plus_back = wavemodel.TransformFromSpatialDomainWithG( w );
w_error = error2(wavemodel.w_plus,w_plus_back);
fprintf('\tG-transform: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(w_error)))))));

% Check the baroclinic F transform
u = wavemodel.TransformToSpatialDomainWithF( wavemodel.u_plus );
u_plus_back = wavemodel.TransformFromSpatialDomainWithF( u );

% In the wave model we inserted an imaginary number at k=l=0 in order to
% randomize the phase... this doesn't have a fourier transform, so we need
% to get rid of it to do a proper comparision.
u_model_fixed = wavemodel.u_plus;
u_model_fixed(1,1,:) = real(u_model_fixed(1,1,:));
u_error = error2(u_model_fixed,u_plus_back);
fprintf('\tF-transform: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(u_error)))))));

uu = wavemodel.TransformToSpatialDomainWithBarotropicFMode(B0);
B0uu = wavemodel.TransformFromSpatialDomainWithBarotropicFMode(uu);
uu_error = error2(B0uu,B0);
fprintf('\tBarotropic transform: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(uu_error)))))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Decomposition test
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n********** Decomposition tests **********\n');

if shouldUseArbitraryStratificationModel == 0
    newmodel = InternalWaveModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);
else
    newmodel = InternalWaveModelArbitraryStratification([Lx, Ly, Lz], [Nx, Ny, Nz], rho, z, latitude, 'method','spectral','cacheFile','evp-cache.mat'); %, 'method','spectral'
end

error2 = @(u,u_unit) max(max(max( abs((u(abs(u_unit)>1e-15)-u_unit(abs(u_unit)>1e-15)))./abs(u_unit(abs(u_unit)>1e-15)) )));

t = 360;
for i=1:8
    if i <= 4 % First we walk through the four types of solutions in isolation
        mask = zeros(4,1);
        mask(i)=1;
    elseif i <= 8  % now walk through triplets
        mask = ones(4,1);
        mask(i-4) = 0;
    else
        mask = ones(4,1);
    end
    wavemodel.A0 = mask(1)*A0;
    wavemodel.GenerateWavePhases(mask(2)*Ap,mask(2)*Am);
    wavemodel.GenerateGeostrophicCurrents(mask(3)*B0,mask(4)*B);
    [u,v,w,eta] = wavemodel.VariableFieldsAtTime(t,'u','v','w','zeta');
    newmodel.InitializeWithHorizontalVelocityAndIsopycnalDisplacementFields(t,u,v,eta);
    
    fprintf('\nmask %d\n',i);
    
    spectralEnergy = 0;
    if ~isempty(error2(newmodel.A0,wavemodel.A0))
        fprintf('The A0 amplitude matches to 1 part in 10^%d\n', round((log10( error2(newmodel.A0,wavemodel.A0) ))));
        % 1 == newmodel.Apm_HKE_factor + newmodel.Apm_VKE_factor + newmodel.Apm_PE_factor
        spectralEnergy = spectralEnergy + abs(newmodel.A0)^2*newmodel.Lz/2;
    end
    if ~isempty(error2(newmodel.Amp_plus,wavemodel.Amp_plus))
        fprintf('The A_plus amplitude matches to 1 part in 10^%d\n', round((log10( error2(newmodel.Amp_plus,wavemodel.Amp_plus) ))));
        % 1 == newmodel.Apm_HKE_factor + newmodel.Apm_VKE_factor + newmodel.Apm_PE_factor
        spectralEnergy = spectralEnergy + sum(sum(sum( newmodel.Amp_plus.*conj(newmodel.Amp_plus) )));
    end
    if ~isempty(error2(newmodel.Amp_minus,wavemodel.Amp_minus))
        fprintf('The A_minus amplitude matches to 1 part in 10^%d\n', round((log10( error2(newmodel.Amp_minus,wavemodel.Amp_minus) ))));
        spectralEnergy = spectralEnergy + sum(sum(sum( newmodel.Amp_minus.*conj(newmodel.Amp_minus) )));
    end
    if ~isempty(error2(newmodel.B0,wavemodel.B0))
        fprintf('The B0 amplitude matches to 1 part in 10^%d\n', round((log10( error2(newmodel.B0,wavemodel.B0) ))));
        spectralEnergy = spectralEnergy + sum(sum(sum( newmodel.B0_HKE_factor .* newmodel.B0.*conj(newmodel.B0) )));
    end
    if ~isempty(error2(newmodel.B,wavemodel.B))
        fprintf('The B amplitude matches to 1 part in 10^%d\n', round((log10( error2(newmodel.B,wavemodel.B) ))));
        spectralEnergy = spectralEnergy + sum(sum(sum( (newmodel.B_HKE_factor+newmodel.B_PE_factor) .* newmodel.B.*conj(newmodel.B) )));
    end
    
    integratedEnergy = trapz(wavemodel.z,mean(mean( u.^2 + v.^2 + w.^2 + shiftdim(wavemodel.N2,-2).*eta.*eta, 1 ),2 ) )/2;
    fprintf('total integrated energy: %f m^3/s\n', integratedEnergy);
    fprintf('total spectral energy: %f m^3/s\n', spectralEnergy);
end

return

% error = @(u,u_unit) max(max(max(abs((u-u_unit)./max(abs(u_unit),1e-100)))));
% A_plus_error = error(newmodel.Amp_plus,wavemodel.Amp_plus);

error2 = @(u,u_unit) abs((u-u_unit))./abs(u_unit);
totalError = error2(newmodel.Amp_plus,wavemodel.Amp_plus);
totalError(abs(wavemodel.Amp_plus)<1e-15) = 0;
A_plus_error = max(max(max(totalError)));

fprintf('The A_plus amplitude matches to 1 part in 10^%d\n', round((log10(A_plus_error))));

totalError = error2(newmodel.Amp_minus,wavemodel.Amp_minus);
totalError(abs(wavemodel.Amp_minus)<1e-15) = 0;
A_minus_error = max(max(max(totalError)));
fprintf('The A_minus amplitude matches to 1 part in 10^%d\n', round((log10(A_minus_error))));

totalError = error2(newmodel.B0,wavemodel.B0);
totalError(abs(wavemodel.B0)<1e-15) = 0;
B0_error = max(max(max(totalError)));
fprintf('The B0 amplitude matches to 1 part in 10^%d\n', round((log10(B0_error))));

totalError = error2(newmodel.B,wavemodel.B);
totalError(abs(wavemodel.B)<1e-15) = 0;
B_error = max(max(max(totalError)));
fprintf('The B amplitude matches to 1 part in 10^%d\n', round((log10(B_error))));

fprintf('Total B amplitude is %f\n', sum(sum(sum(abs(wavemodel.B).^2))) / sum(sum(sum(abs(wavemodel.Amp_minus).^2))))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Depth integrated energy comparison
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


integratedEnergy = trapz(wavemodel.z,mean(mean( u.^2 + v.^2 + w.^2 + N0*N0*eta.*eta, 1 ),2 ) )/2;
fprintf('total integrated energy: %f m^3/s\n', integratedEnergy);

P2_pm = newmodel.Amp_plus.*conj(newmodel.Amp_plus) + newmodel.Amp_minus.*conj(newmodel.Amp_minus);
spectralEnergy = sum(sum(sum(newmodel.Amp_plus.*conj(newmodel.Amp_plus) + newmodel.Amp_minus.*conj(newmodel.Amp_minus))));
fprintf('total spectral energy: %f m^3/s\n', spectralEnergy);

omega = wavemodel.Omega;
f = wavemodel.f0;
N = N0;

HKE = sum(sum(sum(  P2_pm .* (1 + f*f./(omega.*omega)) .* (N*N - omega.*omega) / (2 * (N*N - f*f) )  )));
VKE = sum(sum(sum(  P2_pm .* (omega.*omega - f*f) / (2 * (N*N - f*f) )  )));
PE = sum(sum(sum(  P2_pm .* N*N .* (omega.*omega - f*f) ./ (2 * (N*N - f*f) * omega.*omega )  )));

fprintf('total spectral energy (HKE + VKE + PE) = E: (%f + %f + %f) = %f m^3/s\n', HKE, VKE, PE, HKE+VKE+PE);

return




return

u_in = randn(size(u));
v_in = randn(size(v));

eta_in = randn(size(eta));
eta_in(:,:,1) = 0;
eta_in(:,:,end) = 0;
eta_in = eta_in - mean(mean(eta_in,1),2);

newmodel.InitializeWithHorizontalVelocityAndIsopycnalDisplacementFields(0,u_in,v_in,eta_in);
[u_out,v_out,eta_out] = wavemodel.VariableFieldsAtTime(t,'u','v','zeta');

totalError = error2(u_out,u_in);
totalError(abs(u_in)<1e-15) = 0;
u_error = max(max(max(totalError)));
fprintf('The u amplitude matches to 1 part in 10^%d\n', round((log10(u_error))));

totalError = error2(v_out,v_in);
totalError(abs(v_in)<1e-15) = 0;
v_error = max(max(max(totalError)));
fprintf('The v amplitude matches to 1 part in 10^%d\n', round((log10(v_error))));

totalError = error2(eta_out,eta_in);
totalError(abs(eta_in)<1e-15) = 0;
eta_error = max(max(max(totalError)));
fprintf('The eta amplitude matches to 1 part in 10^%d\n', round((log10(eta_error))));

return;

ubar = wavemodel.TransformFromSpatialDomainWithF( u )./wavemodel.F;
vbar = wavemodel.TransformFromSpatialDomainWithF( v )./wavemodel.F;
etabar = wavemodel.TransformFromSpatialDomainWithG( eta )./wavemodel.G;

k = wavemodel.K;
l = wavemodel.L;
K = sqrt(k.*k + l.*l);
omega = (wavemodel.Omega);
h = wavemodel.h;
g = wavemodel.g;
f0 = wavemodel.f0;

delta = sqrt(h).*(k.*ubar + l.*vbar)./K;
zeta = sqrt(h).*(k.*vbar - l.*ubar)./K;

P_plus = exp(-sqrt(-1)*omega*t).*(-g*K.*sqrt(h).*etabar./omega + delta - sqrt(-1)*zeta*f0./omega)/2;
P_minus = exp(-sqrt(-1)*omega*t).*(g*K.*sqrt(h).*etabar./omega + delta + sqrt(-1)*zeta*f0./omega)/2;
B = (etabar*f0 - sqrt(-1)*zeta.*K.*sqrt(h))*f0./(omega.*omega);

% inertial must be solved for separately.
P_plus(1,1,:) = exp(-sqrt(-1)*f0*t)*(ubar(1,1,:) - sqrt(-1)*vbar(1,1,:)).*sqrt(h(1,1,:))/2;
P_minus(1,1,:) = conj(P_plus(1,1,:));
B(1,1,:) = 0;