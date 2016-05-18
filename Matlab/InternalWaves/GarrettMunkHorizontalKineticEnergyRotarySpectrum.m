%GarrettMunkHorizontalKineticEnergyRotarySpectrum
%
% Given a range of frequencies (in radians) and a buoyancy frequency N,
% this returns the horizontal kinetic energy spectrum.

function [S_gm] = GarrettMunkHorizontalKineticEnergyRotarySpectrum( omega, latitude, z, rho, j_star, drifterDepth, shouldAlias )

b = 1.3e3; % thermocline exponential scale, meters
N0 = 5.2e-3; % reference buoyancy frequency, 1/seconds
E0 = 6.3e-5; % non-dimensional energy parameter
f0 = abs(corfreq(latitude)/3600);
GM_scale = 1.0; % 

Gamma = abs(HKEVerticalStructureFunction( z, rho, latitude, j_star ));
if z(2)-z(1) < 1
    depthIndex = find(z<drifterDepth,1,'first');
else
    depthIndex = find(z<drifterDepth,1,'last');
end

fprintf('Gamma at %f is %f', z(depthIndex), Gamma(depthIndex));

% Regular
if (shouldAlias ~= 1 )
    prefactor = Gamma(depthIndex)*abs(real(GM_scale*E0*b*b*b*N0*N0*(2/pi)*(f0./omega).*(1./sqrt(omega.*omega-f0*f0))));
    S_gm = real(prefactor.*((f0 - omega).^2)./(omega.*omega));
else 
    % Aliased
    N=length(omega);
    M=5;
    sampleInterval = omega(2)-omega(1);
    On = sampleInterval*((-M*N/2):1:(M*N/2)-1)';
    prefactor = Gamma(depthIndex)*abs(real(GM_scale*E0*b*b*b*N0*N0*(2/pi)*(f0./On).*(1./sqrt(On.*On-f0*f0))));
    S_gm = real(prefactor.*((f0 - On).^2)./(On.*On));
    S_gm(find(isnan(S_gm)))=0;

    S_gm = reshape(S_gm, [N M]);
    S_gm = sum(S_gm,2);
end
