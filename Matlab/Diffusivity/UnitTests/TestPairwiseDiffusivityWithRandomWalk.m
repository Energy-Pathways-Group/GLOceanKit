Nx = 10;
Ny = 10;
reps = Nx*Ny;
kappa = .1;

t = (0:1000)';
deltaT = t(2)-t(1);
N = length(t);
sigma = sqrt(2*kappa/deltaT);

dX = sigma*randn(N,reps);
dY = sigma*randn(N,reps);

rmsDisplacement = sqrt(kappa*max(t))

% If I don't put these in a grid, then I get a length scale dependent
% diffusivity. Why? Because I'm grabbing and sorting, biasing the sample.
dL = 10; % sets the artificial grid space

[x0, y0] = ndgrid( dL*(1:Nx), dL*(1:Ny) );
x0 = reshape(x0, 1, reps);
y0 = reshape(y0, 1, reps);

x = deltaT*cumsum(dX,1) + repmat(x0,N,1);
y = deltaT*cumsum(dY,1) + repmat(y0,N,1);

figure
[r2, kappa_r_slope, ~] = PairwiseRelativeDiffusivity( t, x, y, 'slope');

theBins = dL*(1:Nx) + dL/2;
histogramWithErrorbars(sqrt(r2),kappa_r_slope,theBins);
ylim([0 1.2*max(yMean+yStdErr)])

hold on
[r2, kappa_r_powspec, ~] = PairwiseRelativeDiffusivity( t, x, y, 'powspec');
histogramWithErrorbars(sqrt(r2),kappa_r_powspec,theBins);

[r2, kappa_r_endpoint, ~] = PairwiseRelativeDiffusivity( t, x, y, 'endpoint');
histogramWithErrorbars(sqrt(r2),kappa_r_endpoint,theBins);

[r2, kappa_r_yv, ~] = PairwiseRelativeDiffusivity( t, x, y, 'yv');
histogramWithErrorbars(sqrt(r2),kappa_r_yv,theBins);

legend('slope','powspec','endpoint','yv')