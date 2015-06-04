function [x_com, y_com, q, r] = CenterOfMass( x, y )
%% CenterOfMass
%
%	Jeffrey J. Early, 2014
%
% Compute the center of mass of a number of particles, (x_com, y_com).
% The function will optionally return the particles in the center of mass reference frame, (q,r)
%
% [x_com, y_com] = CenterOfMass( x, y )
%
% [x_com, y_com, q, r] = CenterOfMass( x, y )

nTime = size(x,1);
nDrifters = size(x,2);

% Compute the center of mass (com)
x_com = mean(x,2);
y_com = mean(y,2);

% Convert drifters into center of mass coordinates (q,r)
q = x - repmat(x_com,[1 nDrifters]);
r = y - repmat(y_com,[1 nDrifters]);