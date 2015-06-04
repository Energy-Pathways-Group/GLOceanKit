function [X2, Y2] = MeanSquareSeparation( x, y )
%% MeanSquareSeparation
%
%	Jeffrey J. Early, 2014
%
% Compute the mean-square separation of the particles by brute force.
% (x,y) are matrices where each row in a time point and each column is a particle.
%
% Note that when computing diffusivity *each* direction grows at a rate of 2*kappa.
% So if you add the two together, M2=X2+Y2, you want 1/4 of the growth rate.

nTime = size(x,1);
nDrifters = size(x,2);

X2 = zeros(nTime,1);
Y2 = zeros(nTime,1);

for iDrifter=1:nDrifters
	for jDrifter=1:nDrifters
		if (iDrifter ~= jDrifter)
			X2 = X2 + (x(:,iDrifter)-x(:,jDrifter)).^2;
			Y2 = Y2 + (y(:,iDrifter)-y(:,jDrifter)).^2;
		end
	end
end

% This is LaCasce (2008) equation 3.
% He uses 2*N*(N-1) as a divisor. The N-1 factor is used in statistics, according to his
% footnote. Instead, we just use N, to be consistent with our physical definition.
% Obviously the difference is negligible for large N.
X2 = X2 ./ (2*nDrifters*(nDrifters-0));
Y2 = Y2 ./ (2*nDrifters*(nDrifters-0));