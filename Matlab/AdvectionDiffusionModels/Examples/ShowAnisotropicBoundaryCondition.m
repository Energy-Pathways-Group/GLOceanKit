reps = 1000;
nDims = 2;

t = (0:5:9000)';
deltaT = 1;
N = length(t);

% Notation is that f is the flux, which gets dt, and g gets dW
% first column, x, second column u
%
% dx = (\bar{u} + u)*dt
% du = (-(u/T)+sigma_x)*dt + sigma*dW_t

tau = 12; % time scale of the OU
L = 100; % box width
sigma0 = (.1)^2; % velocity variance
sigma = @(x) sigma0 * (-4*x.*x/L^2 + 4*x/L); % velocity variance as function of x
halfsigma_x = @(x) sigma0 * (-4*x/L^2 + 2/L); % half the derivative of the velocity variance

% the acceleration term. Note that equation 16 in Berloff is incorrectly
% copied from their appendix B... multiplying the 2nd term by 2 is
% incorrect as confirmed by 
a = @(x,u) halfsigma_x(x) + halfsigma_x(x).*u.*u./sigma(x);

f = @(t,x) cat(2,x(:,2),-x(:,2)/tau + a(x(:,1),x(:,2)) );
g = @(t,x) cat(2,zeros(size(x(:,1))), sqrt(2*sigma(x(:,1))/tau) );

x0 = cat(2,linspace(.1,L-.1,reps).',zeros(reps,1));

integrator = IntegratorEulerMaruyama( f, g, x0, 0.1*deltaT );
pn = integrator.IntegrateAlongDimension(t);
x = squeeze(pn(:,1,:)).';
u = squeeze(pn(:,2,:)).';

D2 = x(end,:).^2 ;
kappa_out = mean(D2)/(2*t(end))

figure
subplot(2,1,1)
plot(t,x)
ylim([0 L])
subplot(2,1,2)
histogram(x(end,:))