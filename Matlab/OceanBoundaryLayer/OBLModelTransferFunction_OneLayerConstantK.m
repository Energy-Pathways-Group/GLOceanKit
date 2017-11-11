% This implements model 1b from S. Elipot and S. T. Gille: Ekman layers in the Southern Ocean.
% www.ocean-sci.net/5/115/2009
%
% Time is in days. Depth in meters. Latitude in degrees. h is the layer depth. K0 has units of m^2/s
% The returned transfer function has units of m^2 s / kg
function [H] = OBLModelTransferFunction_OneLayerConstantK( time, z, latitude, h, K0 );

f = (2 * 7.2921e-5)*sind(latitude);
rho_water = 1025; % units of kg/m^3

% Lets go to the frequency domain
% nu has units of 1/s
nu = FFTFrequenciesFromTimeSeries( time*86400 );

% Make a an appropriate coordinate grid.
% The resulting matrices are length(z) x length(nu).
[NU, Z] = meshgrid( nu, z);

% Scale depth for each frequency
delta_1 = sqrt(2*K0./(2*pi*NU+f)); % units of m

% Transfer function, units of m^2 s / kg
H = exp( -sqrt(-1)*pi/4 ) ./ ( rho_water * sqrt( (2*pi*NU + f)*K0 ) );
H = H .* sinh( (1+sqrt(-1))*(h-Z)./delta_1 ) ./ cosh( (1+sqrt(-1))*h./delta_1 );

