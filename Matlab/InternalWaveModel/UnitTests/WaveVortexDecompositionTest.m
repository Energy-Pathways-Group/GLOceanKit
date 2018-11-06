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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if shouldUseArbitraryStratificationModel == 0
    wavemodel = InternalWaveModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);
else
    if ~exist('wavemodel','var')
        wavemodel = InternalWaveModelArbitraryStratification([Lx, Ly, Lz], [Nx, Ny, Nz], rho, z, latitude);
    end
end

wavemodel.InitializeWithGMSpectrum(1.0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Forward/back transformation tests
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error = @(u,u_unit) max( [max(max(max(abs(u-u_unit)/max( [max(max(max( abs(u) ))), 1e-15] )))), 1e-15]);
error2 = @(u,u_unit) abs((u-u_unit))./(max(max(max(abs(u_unit)))));

w = wavemodel.TransformToSpatialDomainWithG( wavemodel.w_plus );
w_plus_back = wavemodel.TransformFromSpatialDomainWithG( w );

w_error = error2(wavemodel.w_plus,w_plus_back);
fprintf('The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(w_error)))))));

u = wavemodel.TransformToSpatialDomainWithF( wavemodel.u_plus );
u_plus_back = wavemodel.TransformFromSpatialDomainWithF( u );

% In the wave model we inserted an imaginary number at k=l=0 in order to
% randomize the phase... this doesn't have a fourier transform, so we need
% to get rid of it to do a proper comparision.
u_model_fixed = wavemodel.u_plus;
u_model_fixed(1,1,:) = real(u_model_fixed(1,1,:));
u_error = error2(u_model_fixed,u_plus_back);
fprintf('The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(u_error)))))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Decomposition test
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = 360;
[u,v,w,eta] = wavemodel.VariableFieldsAtTime(t,'u','v','w','zeta');

if shouldUseArbitraryStratificationModel == 0
    newmodel = InternalWaveModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);
else
    if ~exist('newmodel','var')
        newmodel = InternalWaveModelArbitraryStratification([Lx, Ly, Lz], [Nx, Ny, Nz], rho, z, latitude);
    end
end

newmodel.InitializeWithHorizontalVelocityAndIsopycnalDisplacementFields(t,u,v,eta);


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