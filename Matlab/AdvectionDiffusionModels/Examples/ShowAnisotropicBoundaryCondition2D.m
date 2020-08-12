Nx = 10;
t = (0:5:900)'; % total integration time/dimension
deltaT = .01; % actual time step

% Notation is that f is the flux, which gets dt, and g gets dW
% first column, x, second column y, third column u, fourth column v
%
% dx = (\bar{u} + u)*dt
% du = (-(u/T)+sigma_x)*dt + sigma*dW_t

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Parameters
%
L = 100; % box width

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OU and associated parameters
%
tau = 20; % time scale of the OU
sigma0 = (.1)^2; % velocity variance
sigma1 = @(x) sigma0 * sin(pi*x/L); % velocity variance as function of x
halfsigma1_x = @(x) sigma0 * pi * cos(pi*x/L)/(2*L); % half the derivative of the velocity variance
sigma2 = @(y) sigma0 * sin(pi*y/L); % velocity variance as function of x
halfsigma2_y = @(y) sigma0 * pi * cos(pi*y/L)/(2*L); % half the derivative of the velocity variance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Background velocity and associated parameters
%
U0 = 0.05;
U = @(x,y) -U0*pi*sin(pi*x/L).*cos(pi*y/L);
V = @(x,y) U0*pi*cos(pi*x/L).*sin(pi*y/L);


% the acceleration term. Note that equation 16 in Berloff is incorrectly
% copied from their appendix B... multiplying the 2nd term by 2 is
% incorrect, confirmed numerically.
a1 = @(x,y,u) halfsigma1_x(x) + halfsigma1_x(x).*(u+U(x,y)).*u./sigma1(x);
a2 = @(x,y,v) halfsigma2_y(y) + halfsigma2_y(y).*(v+V(x,y)).*v./sigma2(y);

% A little messy, but x=x(:,1),y=x(:,2),u=x(:,3),v=x(:,4)
f = @(t,x) cat(2,U(x(:,1),x(:,2))+x(:,3),V(x(:,1),x(:,2))+x(:,4),-x(:,3)/tau + a1(x(:,1),x(:,2),x(:,3)),-x(:,4)/tau + a2(x(:,1),x(:,2),x(:,4)) );
g = @(t,x) cat(2,zeros(size(x(:,1))),zeros(size(x(:,2))), sqrt(2*sigma1(x(:,1))/tau), sqrt(2*sigma2(x(:,2))/tau) );

% Now place particles evenly spaced
x0 = linspace(.1,L-.1,Nx).';
y0 = linspace(.1,L-.1,Nx).';
[X,Y] = ndgrid(x0,y0);
x0 = X(:);
y0 = Y(:);

x0 = cat(2,x0,y0,zeros(length(x0),1),zeros(length(x0),1));

integrator = IntegratorEulerMaruyama( f, g, x0, deltaT );
pn = integrator.IntegrateAlongDimension(t);
x = squeeze(pn(:,1,:)).';
y = squeeze(pn(:,2,:)).';

figure
plot(x,y)
