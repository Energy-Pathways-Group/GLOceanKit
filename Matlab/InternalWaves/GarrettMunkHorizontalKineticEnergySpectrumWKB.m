%GarrettMunkHorizontalKineticEnergyRotarySpectrum
%
% Given a range of frequencies (in radians) and a buoyancy frequency N,
% this returns the horizontal kinetic energy spectrum.

function [S_gm] = GarrettMunkHorizontalKineticEnergySpectrumWKB( omega, latitude, N_depth, shouldAlias )

b = 1.3e3; % thermocline exponential scale, meters
N0 = 5.2e-3; % reference buoyancy frequency, radians/seconds
E0 = 6.3e-5; % non-dimensional energy parameter
f0 = abs(corfreq(latitude)/3600);
GM_scale = 1.0; % 

% B = (f0./omega).*(1./sqrt(omega.*omega-f0*f0));
% prefactor = GM_scale*(b*b*N*N0*E0/pi)*B.*((N*N - omega.*omega)/(N*N));
% snn_gm = real(prefactor.*((omega+f0).^2)./(omega.*omega));
% spp_gm = real(prefactor.*((omega-f0).^2)./(omega.*omega));
% snn_gm(find(omega<f0))=1e-7;
% spp_gm(find(omega<f0))=1e-7;

% Regular

% I don't like this extra factor because it makes the sum not integrable...
% and I don't know where it comes from.
lienFactor = ((N_depth*N_depth - omega.*omega)/(N_depth*N_depth));

if (shouldAlias ~= 1 )
    %prefactor = abs(real(GM_scale*(b*b*N_depth*N0*E0/pi)*(f0./omega).*(1./sqrt(omega.*omega-f0*f0)).*lienFactor));
    prefactor = abs(real(GM_scale*(b*b*N_depth*N0*E0/pi)*(f0./omega).*(1./sqrt(omega.*omega-f0*f0))));
    S_gm = real(prefactor.*((f0*f0 + omega.*omega))./(omega.*omega));
else 
    % % Aliased
    N=length(omega);
    M=51;
    sampleInterval = omega(2)-omega(1);
    On = sampleInterval*((-M*N/2):1:(M*N/2)-1)';
    %prefactor = abs(real(GM_scale*(b*b*N_depth*N0*E0/pi)*(f0./On).*(1./sqrt(On.*On-f0*f0)).*lienFactor));
    prefactor = abs(real(GM_scale*(b*b*N_depth*N0*E0/pi)*(f0./On).*(1./sqrt(On.*On-f0*f0))));
    S_gm = real(prefactor.*((f0 - On).^2)./(On.*On));
    S_gm(find(isnan(S_gm)))=0;

    S_gm = reshape(S_gm, [N M]);
    S_gm = sum(S_gm,2);
end

% 
% % Blurred
% MF = floor(N/2)+1;
% SpecPermute = [S_gm(MF:N); S_gm(1:MF-1)]; % Permute Spectrum for ifft
% Autocov = ifft(SpecPermute); % Autocovariance of Spectrum
% 
% TriangleKernel = (1-([0:MF-1 ceil(N/2)-1:-1:1]/N))'; % Triange kernel (works for odd and even N)
% S_gm = fftshift(real(fft(TriangleKernel.*Autocov))); % Blurred Spectrum