Nx = 10;
Ny = 10;
t = (0:60:86400)';
deltaT = t(2)-t(1);

xoffset = -1*dL*Nx/2;
yoffset = -1*dL*Ny/2;

dL = 100; % sets the artificial grid space
[x0, y0] = ndgrid( dL*(1:Nx)+xoffset, dL*(1:Ny)+yoffset );
x0 = reshape(x0, [], 1);
y0 = reshape(y0, [], 1);

kappa = 0.00;
sigma = 1e-5;


% If I don't put these in a grid, then I get a length scale dependent
% diffusivity. Why? Because I'm grabbing and sorting, biasing the sample.
[x, y] = ParticlePathInStrainVorticityField( x0, y0, t, 0, sigma, 0, 0, 0, kappa );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%
subplot(2,1,1)
%%%%%%%%%%%%%%
[r2, kappa_r_slope] = PairwiseRelativeDiffusivity( t, x, y, 'slope');

theBins = dL*(1:Nx) + dL/2;
[xMean, yMean, yStdErr] = histogramWithErrorbars(sqrt(r2),kappa_r_slope,theBins);
ylim([0 1.2*max(yMean+yStdErr)])

hold on
[r2, kappa_r_powspec, ~] = PairwiseRelativeDiffusivity( t, x, y, 'powspec');
histogramWithErrorbars(sqrt(r2),kappa_r_powspec,theBins);

[r2, kappa_r_endpoint, ~] = PairwiseRelativeDiffusivity( t, x, y, 'endpoint');
histogramWithErrorbars(sqrt(r2),kappa_r_endpoint,theBins);

[r2, kappa_r_yv, ~] = PairwiseRelativeDiffusivity( t, x, y, 'yv');
histogramWithErrorbars(sqrt(r2),kappa_r_yv,theBins);

title('Relative diffusivity (strain field)')
legend('slope','powspec','endpoint','yv')
xlog, ylog
xlim([min(xMean) max(xMean)])

%%%%%%%%%%%%%%
subplot(2,1,2)
%%%%%%%%%%%%%%
[r2, kappa_a_slope] = SingleParticleDiffusivity( t, x, y, 'slope');

theBins = dL*(1:Nx) + dL/2;
[xMean, yMean, yStdErr] = histogramWithErrorbars(sqrt(r2),kappa_a_slope,theBins);
ylim([0 1.2*max(yMean+yStdErr)])

hold on
[r2, kappa_a_powspec] = SingleParticleDiffusivity( t, x, y, 'powspec');
[xMean, yMean, yStdErr] = histogramWithErrorbars(sqrt(r2),kappa_a_powspec,theBins);

[r2, kappa_a_endpoint] = SingleParticleDiffusivity( t, x, y, 'endpoint');
histogramWithErrorbars(sqrt(r2),kappa_a_endpoint,theBins);

[r2, kappa_a_yv] = SingleParticleDiffusivity( t, x, y, 'yv');
histogramWithErrorbars(sqrt(r2),kappa_a_yv,theBins);

title('Absolute diffusivity (strain field)')
legend('slope','powspec','endpoint','yv')
xlog, ylog
xlim([min(xMean) max(xMean)])

% Compute the slope of the powspec data
[p,S,mu]=polyfit(log(xMean),log(yMean),1);
slope = p(1)/mu(2);
intercept = p(2)-p(1)*mu(1)/mu(2);
fprintf('Actual: %.2g x^(%.2f)\n', exp(intercept),slope);

% Using my notes "frequency-spectrum-from-strain-field.pdf"
expectedSlope = sigma*tanh(sigma*t(end)/4)/16;
expectedPower = 2;
fprintf('Expected: %.2g x^(%.2f)\n', expectedSlope,expectedPower);