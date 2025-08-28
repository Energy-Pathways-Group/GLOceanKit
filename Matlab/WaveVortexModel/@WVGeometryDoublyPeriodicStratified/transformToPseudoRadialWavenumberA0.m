function [varargout] = transformToPseudoRadialWavenumberA0(self,varargin)     
% transforms in the from (j,kRadial) to kPseudoRadial
%
% Sums all the variance/energy in radial bins `kPseudoRadial`.
%
% - Topic: Operations — Transformations
% - Declaration: [varargout] = transformToRadialWavenumber(varargin) 
% - Parameter varargin: variables with dimensions $$(j,kl)$$
% - Returns varargout: variables with dimensions $$(kRadial)$$ or $$(kRadial,j)$$

% Thi is the final output axis for wavenumber

jWavenumber = 1./sqrt(self.Lr2);
jWavenumber(1) = 0; % barotropic mode is a mean?
% idx = find(wvt.kRadial<jWavenumber(2));
% kPseudoRadial = cat(1,wvt.kRadial(idx),jWavenumber(2:end));
% kPseudoRadial = jWavenumber;
[kj,kr] = ndgrid(jWavenumber,self.kRadial);
Kh = sqrt(kj.^2 + kr.^2);

k = self.kPseudoRadial;
dk = k(2)-k(1);

nK = length(k);

varargout = cell(size(varargin));
spectralMatrixSize = [self.Nj length(self.kRadial)];
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