function fillOutWaveSpectrum(self,maxTimeGap)
% Add external waves to the model to fill out the spectrum
%
% Add free/external waves to fill in the gaps of the gridded
% solution. No gaps will be larger than 2*pi/maxTimeGap, and
% the gaps will be smaller near f.
%
% - Topic: External (non-gridded) modes
if nargin < 2
    maxTimeGap = 86140;
end

% the function dOmegas has an initial gap of dOmegaInitial, and
% asymptotes to maxdOmega. Ramp up puts energy in the
% lowest/most energetic modes.
dOmegaInitial = 0.05*self.f;
maxdOmega = 2*pi/maxTimeGap; % Gaps will appear with observations longer than T = 2*pi/maxdOmega;
Ln = -1/log(1-dOmegaInitial/maxdOmega);
dOmegas = maxdOmega*(1-exp(-(1:100)'/Ln));
gapOmegas = self.f + cumsum(dOmegas);
omegaExt = [];
jExt = [];
Omega = self.Omega;
for iMode = 1:(self.Nj-1)
    omegas = sort(reshape(abs(Omega(:,:,iMode+1)),[],1));

    % First fill in the lower triangle
    indices = find(gapOmegas < omegas(2));
    jExt = cat(1,jExt,iMode*ones(length(indices),1));
    omegaExt = cat(1,omegaExt,gapOmegas(indices));

    % Then fill in overly sized gaps
    diffOmega = diff(omegas);
    gapIndices = find(diffOmega>maxdOmega);
    for i=2:length(gapIndices)
        n = ceil(diffOmega(gapIndices(i))/maxdOmega);
        newOmegas = linspace(omegas(gapIndices(i)),omegas(gapIndices(i)+1),n+1)';
        jExt = cat(1,jExt,iMode*ones(n-1,1));
        omegaExt = cat(1,omegaExt,newOmegas(2:end-1));
    end
end
alphaExt = 2*pi*rand( size(omegaExt) );
phiExt = 2*pi*rand( size(omegaExt) );
UExt = zeros(size(omegaExt));

self.setExternalWavesWithFrequencies(omegaExt,alphaExt,jExt,phiExt,UExt,Normalization.kConstant);

fprintf('Added %d external waves to fill out the GM spectrum.\n', length(omegaExt));
end