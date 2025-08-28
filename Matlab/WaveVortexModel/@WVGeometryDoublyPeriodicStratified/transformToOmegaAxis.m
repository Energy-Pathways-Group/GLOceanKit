function [varargout] = transformToOmegaAxis(self,varargin)     
% transforms in the from (j,kRadial) to omegaAxis
%
% Sums all the variance/energy in radial bins `kPseudoRadial`.
%
% - Topic: Operations — Transformations
% - Declaration: [varargout] = transformToRadialWavenumber(varargin) 
% - Parameter varargin: variables with dimensions $$(j,kl)$$
% - Returns varargout: variables with dimensions $$(kRadial)$$ or $$(kRadial,j)$$

% Thi is the final output axis for wavenumber

wvt = self.wvt;

[omegaN,n] = self.wvt.transformToRadialWavenumber(abs(self.wvt.Omega),ones(size(self.wvt.Omega)));
omegaJK = (omegaN./n);

omega = self.omegaAxis;
dOmega = omega(2)-omega(1);
nOmega = length(omega);

varargout = cell(size(varargin));
spectralMatrixSize = [wvt.Nj length(wvt.kRadial)];
for iVar=1:length(varargin)
    if size(varargin{iVar},2) ~= spectralMatrixSize(2)
        error('The input matrix must be of size [Nj NkRadial]');
    end
    
    varargout{iVar} = zeros([nOmega 1]);
end

totalIndices = false(size(omegaJK));
for iK = 1:1:nOmega
    indicesForK = omega(iK)-dOmega/2 <= omegaJK & omegaJK < omega(iK)+dOmega/2;
    for iVar=1:length(varargin)
        varargout{iVar}(iK) =  sum(varargin{iVar}(indicesForK));
    end
    totalIndices = totalIndices | indicesForK;
end

end