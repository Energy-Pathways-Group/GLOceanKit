function [x, y] = ParticlePathInStrainVorticityField( x0, y0, t, zeta, sigma_n, sigma_s, u_0, v_0, kappa )
%% ParticlePathInStrainVorticityField
%
%	Jeffrey J. Early, 2014
%
% (x0, y0) is the initial position of the particle, in meters.
% This may be a vector for multiple particles.
% t is the time vector, in seconds.
% zeta is the vorticity, in 1/seconds
% sigma_n and sigma_s are the normal and shear strain rate, in 1/seconds
% (u_0, v_0) are the mean velocity, in m/s
% kappa is the diffusivity, in m^2/s

% Force (x0, y0) to be row vectors [1 m]
if (~isrow(x0))
	x0=x0';
end
if (~isrow(y0))
	y0=y0';
end
% Force t to be a column vector [n 1]
if (isrow(t))
	t=t';
end

% Determine the step size for the random walk necessary to produce kappa
deltaT = t(2)-t(1);
N = length(t);
stepsize = sqrt(deltaT*2*kappa);

% Generate the random increments
dX = stepsize*randn(N,length(x0));
dY = stepsize*randn(N,length(y0));

u_0 = u_0*ones(size(x0));
v_0 = v_0*ones(size(y0));

sigma2 = sigma_n*sigma_n + sigma_s*sigma_s;
s2 = sigma2 - zeta*zeta;
s = sqrt(s2);

x = ones( length(t), 1 )*x0;
y = ones( length(t), 1 )*y0;

for n=2:length(t);	
	% row vector, [1 m]
	% the previous position
	xn = x(n-1,:);
	yn = y(n-1,:);
    
    if (sigma2>0 && s2 == 0)
        % The homogenous part of the solution
        x(n,:) = (1+sigma_n*deltaT/2)*xn + ((sigma_s-zeta)*deltaT/2)*yn;
        y(n,:) = ((sigma_s+zeta)*deltaT/2)*xn + (1-sigma_n*deltaT/2)*yn;
        
    elseif (s2 < 0)
        s = sqrt(zeta*zeta-sigma2);
        cos_t = cos(s*deltaT/2);
        sin_t = sin(s*deltaT/2);
        
        % The homogenous part of the solution
		x(n,:) = (cos_t + (sigma_n/s)*sin_t)*xn + ((sigma_s-zeta)/s)*sin_t*yn;
		y(n,:) = ((sigma_s+zeta)/s)*sin_t*xn + (cos_t - (sigma_n/s)*sin_t)*yn;
        
        % The non-homogenous part of the solution
		x(n,:) = x(n,:) + (2/s^2)*( (s*sin_t + sigma_n*(1-cos_t))*u_0 + (sigma_s - zeta)*(1-cos_t)*v_0);
		y(n,:) = y(n,:) + (2/s^2)*( (sigma_s + zeta)*(1-cos_t)*u_0 + (s*sin_t - sigma_n*(1-cos_t))*v_0);
        
    elseif (s2 > 0 )
		% scalar, [1 1]
		cosh_t = cosh(s*deltaT/2);
		sinh_t = sinh(s*deltaT/2);
	
		% The homogenous part of the solution
		x(n,:) = (cosh_t + (sigma_n/s)*sinh_t)*xn + ((sigma_s-zeta)/s)*sinh_t*yn;
		y(n,:) = ((sigma_s+zeta)/s)*sinh_t*xn + (cosh_t - (sigma_n/s)*sinh_t)*yn;

		% The non-homogenous part of the solution
		x(n,:) = x(n,:) + (2/s^2)*( (s*sinh_t + sigma_n*(cosh_t - 1))*u_0 + (sigma_s - zeta)*(cosh_t - 1)*v_0);
		y(n,:) = y(n,:) + (2/s^2)*( (sigma_s + zeta)*(cosh_t - 1)*u_0 + (s*sinh_t - sigma_n*(cosh_t - 1))*v_0);
	else
		x(n,:) = xn;
		y(n,:) = yn;
    end
	
	% The diffusive part of the solution
	x(n,:) = x(n,:) + dX(n,:);
	y(n,:) = y(n,:) + dY(n,:);
end