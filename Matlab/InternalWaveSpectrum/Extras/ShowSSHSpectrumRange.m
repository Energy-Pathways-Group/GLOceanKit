% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 2;
LoadFigureDefaults

% Initialize the GarrettMunkSpectrum class with the profile, but do *not*
% reinitialize it if it already exists. Recomputing the modes can take a
% while.
if ~exist('GM','var')
    GM = GarrettMunkSpectrum('/Users/jearly/Documents/ProjectRepositories/GLOceanKit/Matlab/InternalWaveSpectrum/PrecomputedProfiles/exponential-free-surface');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SSH Spectrum
%
% SWOT error pulled from here: https://swot.jpl.nasa.gov/docs/D-61923_SRD_Rev_A_20160318.pdf
k = exp(linspace(log(2*pi/1e10),log(2*pi/1.5e1),1000));
S = GM.IsopycnalSpectrumAtWavenumbers(0,k);
fprintf('total SSH variance due to IGWs is %.2f cm^2\n', trapz(k,S)*1e4);
% S = S./(2*pi*k);
lambda = k*(1e3)/(2*pi);
S_lambda = 1e4*S*2*pi/1e3;
E = 2 + 0.00125*(lambda).^(-2);

if ~exist('GM_latmix','var')
    GM_latmix = GarrettMunkSpectrum('/Users/jearly/Documents/ProjectRepositories/GLOceanKit/Matlab/InternalWaveSpectrum/PrecomputedProfiles/latmix-site1-free-surface');
end
S_latmix = 2*GM_latmix.IsopycnalSpectrumAtWavenumbers(0,k);
S_latmix_lambda = 1e4*S_latmix*2*pi/1e3;

S1d = zeros(size(S));
S1d_latmix = zeros(size(S));
for i=1:(length(S1d)-1)
   S1d(i) = trapz(log(k(i:end)),S(i:end).* (k(i:end)./sqrt( k(i:end).^2 + k(i)^2)) );
   S1d_latmix(i) = trapz(log(k(i:end)),S_latmix(i:end).* (k(i:end)./sqrt( k(i:end).^2 + k(i)^2)) );
end
fprintf('total SSH variance due to IGWs is %.2f cm^2\n', trapz(k,S1d)*1e4);
S1d_lambda = 1e4*S1d*2*pi/1e3;
S1d_latmix_lambda = 1e4*S1d_latmix*2*pi/1e3;


FigureSize = [50 50 figure_width_2col 400];

fig1 = figure('Units', 'points', 'Position', FigureSize);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];

loglog(lambda,S1d_lambda, 'LineWidth', 1.0*scaleFactor)
hold on
plot(lambda,S_latmix_lambda, 'LineWidth', 1.0*scaleFactor)
plot(lambda,E, 'LineStyle', '--', 'LineWidth', 1.0*scaleFactor, 'Color', 0*[1 1 1])
xlabel('cycle/km', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
ylabel('cm^2/(cycle/km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
title(sprintf('IGW SSH spectrum\n'));
legend(sprintf('1d exp profile, %.2f cm^2 (1 GM)',trapz(k,S)*1e4), sprintf('2d Latmix profile, %.2f cm^2 (2 GM)',trapz(lambda,S1d_latmix_lambda)),'SWOT noise level')
ylim([1 1e4])
xlim([min(1e3/1e6) max(1e3/15e3)])
xticks([1e-3 1e-2 1/16])
labels = cell(3,1); labels{1} = '1000 km';  labels{2} = '100 km';  labels{3} = '15 km';
xticklabels(labels)
set( gca, 'FontSize', figure_axis_tick_size);

% intercept = A^2/lambda_0^(2*alpha)
intercept = 70;
lambda_0 = 1/200;
alpha = 1;
A = sqrt(intercept*lambda_0^(2*alpha));
matern = @(lambda) A^2./(lambda.^2 + lambda_0^2).^alpha;
matern_radial = @(lambda) (pi/2)*A^2* lambda./((lambda.^2 + lambda_0^2).^(alpha+1/2));
hold on, plot(lambda,matern(lambda))
hold on, plot(lambda,matern_radial(lambda))

% print('-depsc2', 'SSHWavenumberSpectralRange.eps')