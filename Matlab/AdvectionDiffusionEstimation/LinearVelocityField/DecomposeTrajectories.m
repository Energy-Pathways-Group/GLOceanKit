function [u_meso,v_meso,u_bg,v_bg,u_sm,v_sm,dmxdt,dmydt] = DecomposeTrajectories(x, y, t, parameterEstimates)
% DecomposeTrajectories Peforms the decomposition of trajectories into
% mesoscale, background and submesoscale, following equations (23)-(25) in
% Oscroft, Sykulski and Early (2020).
%
% Required inputs are,
%   x ? [nT nDrifters] x position in meters of the nDrifters
%   y - [nT nDrifters] y position in meters of the nDrifters
%   t - [nT 1] times in seconds of the observation times
%
% Outputs are,
%   [u_meso,v_meso] - [nT nDrifters] mesoscale velocity in m/s
%   [u_bg,v_bg] - [nT 1] background velocity in m/s
%   [u_sm,v_sm] - [nT nDrifters] submesoscale velocity in m/s
%   [dmxdt,dmydt] - [nT 1] center-of-mass velocity in m/s

u0 = parameterEstimates.u0;
v0 = parameterEstimates.v0;
u1 = parameterEstimates.u1;
v1 = parameterEstimates.v1;
sigma_n = parameterEstimates.sigma_n;
sigma_s = parameterEstimates.sigma_s;
zeta = parameterEstimates.zeta;
delta = parameterEstimates.delta;

u_meso = u0 + u1.*t + 0.5*(sigma_n+delta).*x + 0.5*(sigma_s-zeta).*y;
v_meso = v0 + v1.*t + 0.5*(sigma_s+zeta).*x + 0.5*(delta-sigma_n).*y;

% To compute the background velocity, we need the center-of-mass velocity.
D = FiniteDifferenceMatrix(1,t,1,1,2);
mx = mean(x,2);
my = mean(y,2);
dmxdt = D*mx;
dmydt = D*my;
u_bg = dmxdt - mean(u_meso,2);
v_bg = dmydt - mean(v_meso,2);

% To compute the submesoscale velocity, we need the total velocity as well.
dxdt = D*x;
dydt = D*y;
u_sm = dxdt - u_meso - u_bg;
v_sm = dydt - v_meso - v_bg;

% q = (x-mx);
% r = (y-my);
% dqdt = D*q;
% drdt = D*r;
% dqdt_meso = 0.5*(sigma_n+delta).*q + 0.5*(sigma_s-zeta).*r;
% drdt_meso = 0.5*(sigma_s+zeta).*q + 0.5*(delta-sigma_n).*r;

end

