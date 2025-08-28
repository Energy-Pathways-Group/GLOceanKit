function [varargout] = transformToPseudoRadialWavenumberApm(self,varargin)     
% transforms in the from (j,kRadial) to kPseudoRadial
%
% Sums all the variance/energy in radial bins `kPseudoRadial`.
%
% - Topic: Operations — Transformations
% - Declaration: [varargout] = transformToRadialWavenumber(varargin) 
% - Parameter varargin: variables with dimensions $$(j,kl)$$
% - Returns varargout: variables with dimensions $$(kRadial)$$ or $$(kRadial,j)$$

% Thi is the final output axis for wavenumber

wvt = self.wvt;
if isa(wvt,"WVTransformHydrostatic")
    h_pm = repmat(wvt.h_pm,[1 wvt.Nkl]);
else
    h_pm = wvt.h_pm;
end
[hN,n] = wvt.transformToRadialWavenumber(h_pm,ones(size(h_pm)));
h_kj = (hN./n);

kj = wvt.f ./ sqrt( wvt.g * h_kj);
kj(1) = 0; % barotropic mode is a mean?
kr = repmat( reshape(wvt.kRadial,1,[]),[wvt.Nj 1]);

Kh = sqrt(kj.^2 + kr.^2);

k = self.kPseudoRadial;
dk = k(2)-k(1);

nK = length(k);

varargout = cell(size(varargin));
spectralMatrixSize = [wvt.Nj length(wvt.kRadial)];
for iVar=1:length(varargin)
    if size(varargin{iVar},2) ~= spectralMatrixSize(2)
        error('The input matrix must be of size [Nj NkRadial]');
    end
    
    varargout{iVar} = zeros([nK 1]);
end

totalIndices = false(size(Kh));
for iK = 1:1:nK
    indicesForK = k(iK)-dk/2 <= Kh & Kh < k(iK)+dk/2;
    for iVar=1:length(varargin)
        varargout{iVar}(iK) =  sum(varargin{iVar}(indicesForK));
    end
    totalIndices = totalIndices | indicesForK;
end

end