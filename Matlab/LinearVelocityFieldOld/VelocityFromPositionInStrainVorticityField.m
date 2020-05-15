% Particle position in the center of mass frame. size(x)=[m n] where m is the number
% of time points and n is the number of drifters. size(sigma_n) = [m 1];
function [u, v] = VelocityFromPositionInStrainVorticityField( x, y, zeta, sigma_n, sigma_s )

if numel(sigma_n) == 1 % constant in time
    u = (1/2)*(sigma_n*x + (sigma_s - zeta)*y);
    v = (1/2)*((sigma_s + zeta)*x - sigma_n*y);
elseif size(sigma_n,1) == size(x,1) % time series
    n = size(x,2);
    u = (1/2)*(repmat(sigma_n,[1 n]).*x + (repmat(sigma_s,[1 n]) - repmat(zeta,[1 n])).*y);
    v = (1/2)*((repmat(sigma_s,[1 n]) + repmat(zeta,[1 n])).*x - repmat(sigma_n,[1 n]).*y);
else
   disp('Indice mismatch')
end