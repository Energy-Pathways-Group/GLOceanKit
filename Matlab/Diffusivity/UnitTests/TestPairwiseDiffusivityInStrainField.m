Nx = 10;
Ny = 10;
t = (0:60:86400)';
deltaT = t(2)-t(1);

dL = 100; % sets the artificial grid space
[x0, y0] = ndgrid( dL*(1:Nx), dL*(1:Ny) );
x0 = reshape(x0, 1, reps);
y0 = reshape(y0, 1, reps);

kappa = 0.00;
sigma = 1e-5;


% If I don't put these in a grid, then I get a length scale dependent
% diffusivity. Why? Because I'm grabbing and sorting, biasing the sample.
[x, y] = ParticlePathInStrainVorticityField( x0, y0, t, 0, sigma, 0, 0, 0, kappa );

figure
[r2, kappa_r_slope, ~] = PairwiseRelativeDiffusivity( t, x, y, 'slope');

theBins = dL*(1:Nx) + dL/2;
histogramWithErrorbars(sqrt(r2),kappa_r_slope,theBins);

hold on
[r2, kappa_r_powspec, ~] = PairwiseRelativeDiffusivity( t, x, y, 'powspec');
histogramWithErrorbars(sqrt(r2),kappa_r_powspec,theBins);

[r2, kappa_r_endpoint, ~] = PairwiseRelativeDiffusivity( t, x, y, 'endpoint');
histogramWithErrorbars(sqrt(r2),kappa_r_endpoint,theBins);

[r2, kappa_r_yv, ~] = PairwiseRelativeDiffusivity( t, x, y, 'yv');
histogramWithErrorbars(sqrt(r2),kappa_r_yv,theBins);

legend('slope','powspec','endpoint','yv')