% Mxx0, Myy0, Mxy0 are the initial moments, in units of meters^2.
% t is the time vector, in seconds.
% zeta is the vorticity, 1/seconds
% sigma is the strain, 1/seconds
% theta is the angle of the principal strain axis, a counter-clockwise rotation.
% kappa is the diffusivity, meters^2/second
% Mxx, Myy, Mxy will be returned
function [Mxx, Myy, Mxy] = MomentTensorEvolutionInStrainVorticityField( Mxx0, Myy0, Mxy0, t, zeta, sigma, theta, kappa )

% Force t to be a column vector [n 1]
% if (isrow(t))
% 	t=t';
% end

if sigma==0.0 && zeta == 0.0
	Mxx = 2*kappa*t + Mxx0;
	Myy = 2*kappa*t + Myy0;
	Mxy = 0*t + Mxy0;
elseif zeta==0.0
	cos_t = cos(theta);
	sin_t = sin(theta);

	cos2 = cos_t*cos_t;
	sin2 = sin_t*sin_t;
	cossin = cos_t*sin_t;

	tks = 2*kappa/sigma;
	
	% First take the initial conditions and rotate them to the strain aligned coordinates.
	Mxx0_r = cos2*Mxx0 + sin2*Myy0 + 2*cos_t*sin_t*Mxy0;
	Myy0_r = sin2*Mxx0 + cos2*Myy0 - 2*cos_t*sin_t*Mxy0;
	Mxy0_r = -cossin*Mxx0 + cossin*Myy0 + (cos2 - sin2)*Mxy0;
	
	% Write down the solution in the rotated frame
	Maa = (Mxx0_r + tks)*exp(sigma*t) - tks;
	Mbb = (Myy0_r - tks)*exp(-sigma*t) + tks;
	Mab = Mxy0_r;
	
	% Finally rotate the solution back to the original reference frame.
	Mxx = cos2*Maa + sin2*Mbb - 2*cos_t*sin_t*Mab;
	Myy = sin2*Maa + cos2*Mbb + 2*cos_t*sin_t*Mab;
	Mxy = cossin*Maa - cossin*Mbb + (cos2 - sin2)*Mab;
elseif sigma==0.0
    A = (Mxx0 + Myy0)/2;
    B = -Mxy0;
    C = (Mxx0 - Myy0)/2;
    
    Mxx = A + 2*kappa*t + B*sin(zeta*t) + C*cos(zeta*t);
	Myy = A + 2*kappa*t - B*sin(zeta*t) - C*cos(zeta*t);
	Mxy = -B*cos(zeta*t) + C*sin(zeta*t);
elseif zeta*zeta < sigma*sigma
	cos_t = cos(theta);
	sin_t = sin(theta);

	cos2 = cos_t*cos_t;
	sin2 = sin_t*sin_t;
	cossin = cos_t*sin_t;
	
	s = sqrt( sigma*sigma - zeta*zeta);
	tks = 2*kappa*sigma/(s*s);
	
	% First take the initial conditions and rotate them to the strain aligned coordinates.
	Mxx1 = cos2*Mxx0 + sin2*Myy0 + 2*cos_t*sin_t*Mxy0;
	Myy1 = sin2*Mxx0 + cos2*Myy0 - 2*cos_t*sin_t*Mxy0;
	Mxy1 = -cossin*Mxx0 + cossin*Myy0 + (cos2 - sin2)*Mxy0;
	
	% Compute the coefficients to the solution in the rotate frame.
	A = (1+sigma/s)*Mxx1/2 - (zeta/s)*Mxy1 - (1-sigma/s)*Myy1/2 + tks;
	B = (1-sigma/s)*Mxx1/2 + (zeta/s)*Mxy1 - (1+sigma/s)*Myy1/2 + tks;
	C = -(zeta/s)*Mxx1 + (2*sigma/s)*Mxy1 - (zeta/s)*Myy1;
	
	% Now write down the solution in the rotated frame
	Maa = (A/2)*(1+sigma/s)*exp(s*t) + (B/2)*(1-sigma/s)*exp(-s*t) + (zeta/s)*C/2 - (2*kappa/(s*s))*(zeta*zeta*t + sigma);
	Mbb = -(A/2)*(1-sigma/s)*exp(s*t) - (B/2)*(1+sigma/s)*exp(-s*t) + (zeta/s)*C/2 - (2*kappa/(s*s))*(zeta*zeta*t - sigma);
	Mab = (A/2)*(zeta/s)*exp(s*t) - (B/2)*(zeta/s)*exp(-s*t) + (sigma/s)*C/2 - (2*kappa/(s*s))*(zeta*sigma*t);
	
	% And finally rotate the solution back to the original reference frame.
	Mxx = cos2*Maa + sin2*Mbb - 2*cos_t*sin_t*Mab;
	Myy = sin2*Maa + cos2*Mbb + 2*cos_t*sin_t*Mab;
	Mxy = cossin*Maa - cossin*Mbb + (cos2 - sin2)*Mab;
elseif zeta*zeta == sigma*sigma
    cos_t = cos(theta);
	sin_t = sin(theta);

	cos2 = cos_t*cos_t;
	sin2 = sin_t*sin_t;
	cossin = cos_t*sin_t;
    
    % First take the initial conditions and rotate them to the strain aligned coordinates.
	Mxx1 = cos2*Mxx0 + sin2*Myy0 + 2*cos_t*sin_t*Mxy0;
	Myy1 = sin2*Mxx0 + cos2*Myy0 - 2*cos_t*sin_t*Mxy0;
	Mxy1 = -cossin*Mxx0 + cossin*Myy0 + (cos2 - sin2)*Mxy0;
    
    A = -2*zeta*Mxy1 + sigma*(Mxx1 + Myy1);
    B = Mxx1 - Myy1;
    C = Mxx1 + Myy1;

    % Now write down the solution in the rotated frame
    Maa = kappa*sigma^2*t.^3/3 + kappa*sigma*t.^2 + 2*kappa*t + A*(sigma*t.^2+2*t)/4 + B*(sigma*t+1)/2 + C/2;
    Mbb = kappa*sigma^2*t.^3/3 - kappa*sigma*t.^2 + 2*kappa*t + A*(sigma*t.^2-2*t)/4 + B*(sigma*t-1)/2 + C/2;
    Mab = kappa*sigma*zeta*t.^3/3 + A*(zeta*t.^2/4 - 1/(2*zeta)) + B*zeta*t/2 + sigma*C/(2*zeta);
    
    % And finally rotate the solution back to the original reference frame.
	Mxx = cos2*Maa + sin2*Mbb - 2*cos_t*sin_t*Mab;
	Myy = sin2*Maa + cos2*Mbb + 2*cos_t*sin_t*Mab;
	Mxy = cossin*Maa - cossin*Mbb + (cos2 - sin2)*Mab;
    
elseif zeta*zeta > sigma*sigma
    cos_t = cos(theta);
	sin_t = sin(theta);

	cos2 = cos_t*cos_t;
	sin2 = sin_t*sin_t;
	cossin = cos_t*sin_t;
    
    s = sqrt( zeta*zeta - sigma*sigma);
    
    % First take the initial conditions and rotate them to the strain aligned coordinates.
	Mxx1 = cos2*Mxx0 + sin2*Myy0 + 2*cos_t*sin_t*Mxy0;
	Myy1 = sin2*Mxx0 + cos2*Myy0 - 2*cos_t*sin_t*Mxy0;
	Mxy1 = -cossin*Mxx0 + cossin*Myy0 + (cos2 - sin2)*Mxy0;
    
    A = Mxx1 - Myy1 - 4*kappa*sigma/s^2;
    B = (sigma/s)*Mxx1 - (2*zeta/s)*Mxy1 + (sigma/s)*Myy1;
    C = (zeta^2/s^2)*Mxx1 - (2*zeta*sigma/s^2)*Mxy1 + (zeta^2/s^2)*Myy1;
    
    Maa = (A/2)*(cos(s*t) + (sigma/s)*sin(s*t)) + (B/2)*(sin(s*t)-(sigma/s)*cos(s*t)) + C/2 + (2*kappa/s^2)*(sigma + zeta^2*t);
    Mbb = -(B/2)*(sin(s*t) + (sigma/s)*cos(s*t)) - (A/2)*(cos(s*t)-(sigma/s)*sin(s*t)) + C/2 + (2*kappa/s^2)*(-sigma + zeta^2*t);
    Mab = (A*zeta/(2*s))*sin(s*t) - (B*zeta/(2*s))*cos(s*t) + C*sigma/(2*zeta) + 2*kappa*sigma*zeta*t/s^2;
    
    % And finally rotate the solution back to the original reference frame.
	Mxx = cos2*Maa + sin2*Mbb - 2*cos_t*sin_t*Mab;
	Myy = sin2*Maa + cos2*Mbb + 2*cos_t*sin_t*Mab;
	Mxy = cossin*Maa - cossin*Mbb + (cos2 - sin2)*Mab;
end