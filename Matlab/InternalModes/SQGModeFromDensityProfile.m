%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Given a density profile, this function returns the sqg mode
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [psi, N2] = SQGModeFromDensityProfile( rho, z, latitude, wavelength )

if abs(std(diff(z))/mean(diff(z)))>.000001 
    disp(['Values not evenly spaced!']) 
    return 
end

% We orient the vector so that it is monotonically increasing, and thus
% the bottom is at z(1) and the top at z(end).
didFlip = 0;
if (z(2) - z(1)) < 0
	if (isrow(z))
		z=fliplr(z);
	else
		z=flipud(z);
	end
	if (isrow(rho))
		rho=fliplr(rho);
	else
		rho=flipud(rho);
	end
	didFlip=1;
end

delta = z(2) - z(1);
N = length(z);
g = 9.81;
f0 = corfreq(latitude)/3600;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Construct the the first order differentiation matrix with boundary conditions
% such that F''(endpoints)=0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The operator D is the first order, central difference operator
DiffN = zeros(N,N);
for n=2:N-1
	DiffN(n,n-1) = -1/(2 * delta);
	DiffN(n,n+1) = 1/(2 * delta);
end

% The forwards/backwards difference is used at the end points
DiffN(1,1) = -1/(delta); DiffN(1,2) = 1/(delta);
DiffN(N,N-1) = -1/(delta); DiffN(N,N) = 1/(delta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute the buoyancy
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho_0=mean(rho);

% Convert density to buoyancy. 
Bouyancy=-g*rho/rho_0; 

% N^2 is computed assuming that d^2rho/dz^2=0 at the end points.
N2=DiffN*Bouyancy;

alpha = (f0^2)./N2;

K = 2*pi/wavelength;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Construct the the second order differentiation matrix.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create operator
Diff2 = zeros(N,N);
for n=2:N-1
	Diff2(n,n-1) = alpha(n-1)/(delta^2);
	Diff2(n,n) = -2*alpha(n)/(delta^2) - K*K;
	Diff2(n,n+1) = alpha(n+1)/(delta^2);
end

% bottom of the ocean boundary condition
Diff2(1,1) = - 1 / (delta^1);
Diff2(1,2) = 1 / (delta^1);

% surface of the ocean boundary condition.
Diff2(N,N-1) = -1 / (delta^1);
Diff2(N,N) = 1 / (delta^1);

b = zeros(length(z),1);
b(end)=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Now compute the mode
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

psi = Diff2\b;

% If we flip the variables, flip them back
if (didFlip == 1)
	if (isrow(N2))
		N2=fliplr(N2);
	else
		N2=flipud(N2);
	end
	
	psi = flipdim(psi,1);
end
