
% u_rms = 0.2;
% k_f = wvt.dk*20;
% k_r = wvt.dk*2;
% 
% 
% m = 3/2; % We don't really know what this number is.
% kappa_epsilon = 0.5 * u_rms^2 / ( ((3*m+5)/(2*m+2))*k_r^(-2/3) - k_f^(-2/3) );
% model_viscous = @(k) kappa_epsilon * k_r^(-5/3 - m) * k.^m;
% model_inverse = @(k) kappa_epsilon * k.^(-5/3);
% model_forward = @(k) kappa_epsilon * k_f^(4/3) * k.^(-3);
% model_spectrum = @(k) model_viscous(k) .* (k<k_r) + model_inverse(k) .* (k >= k_r & k<=k_f) + model_forward(k) .* (k>k_f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Setup the model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N0 = 3*2*pi/3600;
L_gm = 1300;
N2 = @(z) N0*N0*exp(2*z/L_gm);
wvt = WVTransformHydrostatic([500e3, 500e3, 4000],[64, 64, 30], N2=N2,latitude=27);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Add initial conditions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Ap,Am,A0] = wvt.geostrophicComponent.randomAmplitudesWithSpectrum(@(k,j) ones(size(k)),@(k,j) ones(size(k)));

Ekj_w = wvt.transformToRadialWavenumber(wvt.A0_TE_factor .* (abs(A0).^2));
figure, plot(wvt.kRadial,sum(Ekj_w,1)), xlog, ylog, xlabel('k'), ylabel('energy')