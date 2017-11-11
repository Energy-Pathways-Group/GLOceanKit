% This implements model 1b from S. Elipot and S. T. Gille: Ekman layers in the Southern Ocean.
% www.ocean-sci.net/5/115/2009
%
% Winds are assumed to be given in meters per second. Time is in days. Depth in meters.
% slab_damp is the slab e-fold damping scale in days. K0 has units of m^2/s
function [t, u_water, v_water] = OBLModel_OneLayerConstantK( time, u_wind, v_wind, z, h, latitude, K0 );

% This returns a complex tau, x+i*y, in N / m^2
[t_stress, tau] = StressFromWindVector( time*86400, u_wind, v_wind);

% Lets go to the frequency domain
% T has units of kg / (m s)
% nu has units of 1/s
[nu, T] = TransformForward( t_stress, tau, 1);

% Make a an appropriate coordinate grid.
% The resulting matrices are length(z) x length(nu).
T = repmat(T.',[length(z) 1]);

% Transfer function, units of m^2 s / kg
H = OBLModelTransferFunction_OneLayerConstantK( time, z, latitude, h, K0 );

U = H .* T; % units of m
[tn, z] = TransformBack( nu, U, 2);

u_water = real(z);
v_water = imag(z);
t = (t_stress/86400);
