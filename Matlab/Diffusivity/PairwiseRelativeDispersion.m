function [D2,r2,r0] = PairwiseRelativeDispersion( t, x, y, theBins )
%PairwiseRelativeDiffusivityFromSlope
%
% Returns the *relative* diffusivity between pairs of drifters as a
% function of initial separation. The results are sorted into bins with
% edges given by "theBins".
%
% Note that the relative diffusivity is twice the usual diffusivity you'd
% get from just computing the rate-of-change of 2nd moment of the mass of
% drifters. So it's double what you'd expect from measuring a tracer.
%
% The algorithm computes the mean-square separation of particle pairs at a
% given time, averages the mean-square, *then* computes the slope of those
% mean separations. The slope is computing assuming an intercept of zero.
%
% This is different than the algorithm used in PairwiseRelativeDiffusivityFromSlopeScatter

nDrifters = size(x,2);

% make sure t is a column vector
t = reshape(t,[],1);

nBins = length(theBins);
D2 = zeros(length(t),nBins); % mean-square distance vs time, for each bin
r0 = zeros(1,nBins); % initial separation
r2 = zeros(1,nBins); % initial separation
nReps = zeros(1,nBins); % number of samples in this bin.

stride = 1;
for iDrifter=1:stride:nDrifters
    for jDrifter = (iDrifter+1):stride:nDrifters        
        q = x(:,iDrifter) - x(:,jDrifter);
        r = y(:,iDrifter) - y(:,jDrifter);
        
        initialSeparation = sqrt(q(1)^2 + r(1)^2);
        iBin = find(initialSeparation > theBins,1,'last');
        
        nReps(iBin) = nReps(iBin) + 1;
        r0(iBin) = r0(iBin) + initialSeparation;
        r2(iBin) = r2(iBin) + initialSeparation^2;
        
        % Now remove the initial conditions.
        q = q-q(1);
        r = r-r(1);
        
        % squared-separation distance as a function of time.
        D2(:,iBin) = D2(:,iBin) + q.^2 + r.^2;
    end
end

D2 = D2./nReps;
r2 = r2./nReps;
r0 = r0./nReps;

