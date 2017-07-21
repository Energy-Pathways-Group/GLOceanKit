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
plot(x,y)
title('Particle positions with initial grid')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now let's go to the center of mass frame.
% Of course, in this case, the center-of-mass isn't moving, so it should be
% nearly zero.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_com = mean(x,2);
y_com = mean(y,2);

% (q,r) are the particle positions in the center of mass frame
q = x - x_com;
r = y - y_com;

% Let's subtract off the initial values (undoing the fake grid we imposed).
% If you don't subtract off the intial conditions, you get very different
% diffusivities, negative even!
q = q - q(1,:);
r = r - r(1,:);

% Now create the second moment matrix [M_qq, M_qr; M_qr, Mrr]
M_qq = mean(q.*q,2);
M_rr = mean(r.*r,2);
M_qr = mean(q.*r,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot what this looks like over time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,2)
plot(t, M_qq + M_rr)

% Now let's fit a line through the second moment
D2 = M_qq + M_rr;
[p,~,mu]=polyfit(t,D2,1);
kappa_fit = (p(1)/mu(2))/4;
intercept = p(2)-p(1)*mu(1)/mu(2);

% Alternatively, we can just use the endpoints
kappa_endpoint = 0.25*(D2(end)-D2(1))/(t(end)-t(1));

hold on
plot(t,4*kappa_fit*t + intercept)
title(sprintf('Second moment (kappa, kappa_{fit}, kappa_{endpoint})=(%.2g,%.2g,%.2g) m^2/s', kappa, kappa_fit, kappa_endpoint))
