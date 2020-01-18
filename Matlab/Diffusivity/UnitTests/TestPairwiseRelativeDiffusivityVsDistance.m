% This generates a grid of particles on a random walk with some specified
% diffusivity.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the grid size, the diffusivity, and the grid spacing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx = 10;
Ny = 10;
kappa = .1;
dL = 10; % sets the artificial grid space

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now generate the random increments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nParticles = Nx*Ny;
t = (0:1000)';
deltaT = t(2)-t(1);
N = length(t);
sigma = sqrt(2*kappa/deltaT);

dX = sigma*randn(N,nParticles);
dY = sigma*randn(N,nParticles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the fake intial conditions (as if the particles are on a grid)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x0, y0] = ndgrid( dL*(1:Nx), dL*(1:Ny) );
x0 = reshape(x0, 1, nParticles);
y0 = reshape(y0, 1, nParticles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sum the increments, and add the initial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = deltaT*cumsum(dX,1) + repmat(x0,N,1);
y = deltaT*cumsum(dY,1) + repmat(y0,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,1)
plot(x,y), axis equal
title('Particle positions with initial grid')

edges = CreateBinEdgesForInitialSeparation(t,x,y);
[r2, kappa_r, kappa_err ] = PairwiseRelativeDiffusivityFromSlope(t, x, y, edges );
[r2_scat, kappa_r_scat, kappa_err_scat ] = PairwiseRelativeDiffusivityFromSlopeScatter(t, x, y, edges );
[r2_pow, kappa_r_pow, kappa_err_pow ] = PairwiseRelativeDiffusivityVsDistance(t, x, y,'powspec',edges);

subplot(1,2,2)
errorbar(sqrt(r2),kappa_r, 2*kappa_err), hold on
errorbar(sqrt(r2_scat),kappa_r_scat, 2*kappa_err_scat)
errorbar(sqrt(r2_pow),kappa_r_pow, 2*kappa_err_pow)
ylim([0 max(kappa_r+2*kappa_err)])
xlabel('distance (m)')
ylabel('kappa (m^2/s)')
legend('slope','scatter-slope','powspec')
