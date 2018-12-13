function [r0,D2] = PairwiseMeanSquareSeparation( t, x, y, edges )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
nDrifters = size(x,2);

% make sure t is a column vector
t = reshape(t,[],1);

nBins = length(edges);
D2 = zeros(length(t),nBins); % mean-square distance vs time, for each bin
r0 = zeros(1,nBins); % initial separation
nReps = zeros(1,nBins); % number of samples in this bin.

nMissing = 0;

stride = 1;
for iDrifter=1:stride:nDrifters
    for jDrifter = (iDrifter+1):stride:nDrifters        
        q = x(:,iDrifter) - x(:,jDrifter);
        r = y(:,iDrifter) - y(:,jDrifter);
        
        initialSeparation = sqrt(q(1)^2 + r(1)^2);
        iBin = find(initialSeparation > edges,1,'last');
        
        if isempty(iBin)
            nMissing = nMissing+1;
            continue;
        end
        
        nReps(iBin) = nReps(iBin) + 1;
        r0(iBin) = r0(iBin) + initialSeparation;
                
        % Now remove the initial conditions.
        q = q-q(1);
        r = r-r(1);
        
        % squared-separation distance as a function of time.
        D2(:,iBin) = D2(:,iBin) + q.^2 + r.^2;
    end
end

if nMissing > 0
    fprintf('%d pairs were not placed into bins.\n',nMissing);
end

D2 = D2./nReps;
r0 = r0./nReps;
end

