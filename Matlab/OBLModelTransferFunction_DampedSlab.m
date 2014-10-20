% This implements the transfer function of a damped slab model.
%
% Time is in days. Depth in meters. Latitude in degrees. 
% slab_damp is the slab e-fold damping scale in days
% The returned transfer function has units of m^2 s / kg
function [H] = OBLModelTransferFunction_DampedSlab( time, latitude, slab_depth, slab_damp )

f = (2 * 7.2921e-5)*sind(latitude);
rho_water = 1025; % units of kg/m^3

if ( isrow(time) )
	time = time.';
end

% Lets go to the frequency domain
% nu has units of 1/s
nu = FFTFrequenciesFromTimeSeries( time*86400 );

% Convert to radians per second
sigma = 2*pi*nu;

% Damping parameter
r = 1/(slab_damp*86400);

% Transfer function, units of m^2 s / kg
H = 1 ./ (rho_water * slab_depth * (r + sqrt(-1) * (f + sigma))) ; % 1/s
%H = 1 ./ (rho_water * slab_depth * ( r + sqrt(-1)*(f + sigma)));

%H = (-r + sqrt(-1)*(f+sigma))/(rho_water*slab_depth);
%H = H ./ (sigma.*sigma - f*f - r*r - 2*sqrt(-1)*r*sigma);