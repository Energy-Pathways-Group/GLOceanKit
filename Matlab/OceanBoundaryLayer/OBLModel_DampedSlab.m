% Winds are assumed to be given in meters per second. Time is in days. Depth in meters.
% slab_damp is the slab e-fold damping scale in days
function [t, u_water, v_water] = OBLModel_DampedSlab( time, u_wind, v_wind, slab_depth, latitude, slab_damp )

if ( isrow(time) )
	time = time.';
end

% This returns a complex tau, x+i*y, in N / m^2
[t_stress, tau] = StressFromWindVector( time*86400, u_wind, v_wind);

if ( isrow(tau) )
	tau = tau.';
end

% Lets go to the frequency domain
% T has units of kg / (m s)
% nu has units of 1/s
[nu, T] = TransformForward( t_stress, tau, 1);

% Transfer function, units of m^2 s / kg
%H = 1; % 1/s
H = OBLModelTransferFunction_DampedSlab( time, latitude, slab_depth, slab_damp);

U = H .* T; % units of m
[tn, z] = TransformBack( nu, U, 1);

u_water = real(z);
v_water = imag(z);
t = (t_stress/86400);
