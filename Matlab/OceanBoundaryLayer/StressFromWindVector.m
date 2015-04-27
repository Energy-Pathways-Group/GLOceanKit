% The wind stress is computed from the wind vector using the empirical coefficients
% found on page 1520 of Yelland, et al., 1998, JPO. Technically this is only valid
% for speeds between 6 and 25 m/s.
function [t, tau] = StressFromWindVector( time, u_wind, v_wind )

uv_wind = u_wind + sqrt(-1)*v_wind;
speed_wind = sqrt(uv_wind .* conj( uv_wind ));

% Dump the nans
%non_nan_indices = find( isnan(speed_wind) ~= 1);
%time = time(non_nan_indices);
%speed_wind = speed_wind(non_nan_indices);
%uv_wind = uv_wind(non_nan_indices);

nan_indices = find( isnan(speed_wind) == 1);
speed_wind(nan_indices) = 0;
uv_wind(nan_indices) = 0;

rho_air = 1.25;
rho_water = 1025;

drag_coefficient = (1e-3)*(0.50 + 0.071*speed_wind);

% Wind stress (not include the speed of the water). Note that uv_wind is not
% a unit vector. Thus we only multiply by one factor of speed_wind.
tau = rho_air * drag_coefficient .* speed_wind .* uv_wind; % N / m^2
t = time;