function [varargout] = transformToRadialWavenumber(self,varargin)     
% transforms in the spectral domain from (j,kl) to (j,kRadial)
%
% Sums all the variance/energy in radial bins `kRadial`.
%
% The following example takes the total energy of the geostrophic part of
% flow, converts it to a one-dimensional spectrum in $$k$$, and then plots
% it with pcolor. The next plot then sums over all wavenumber, and produces
% plots the total energy spectrum as a function of vertical mode $$j$$
% only.
%
% ```matlab
% figure
% tiledlayout('flow')
% Ekj = wvt.transformToRadialWavenumber( wvt.A0_TE_factor .* abs(wvt.A0).^2 );
% nexttile, pcolor(wvt.kRadial,wvt.j,Ekj), shading flat
% nexttile, plot(wvt.j,sum(Ekj,2))
% ```
%
% - Topic: Operations — Transformations
% - Declaration: [varargout] = transformToRadialWavenumber(varargin) 
% - Parameter varargin: variables with dimensions $$(j,kl)$$
% - Returns varargout: variables with dimensions $$(kRadial)$$ or $$(kRadial,j)$$

% Thi is the final output axis for wavenumber
k = self.kRadial;
dk = k(2)-k(1);

nK = length(k);

varargout = cell(size(varargin));
spectralMatrixSize = self.spectralMatrixSize;
for iVar=1:length(varargin)
    if size(varargin{iVar},2) ~= spectralMatrixSize(2)
        error('The input matrix must be of size [Nj Nkl]');
    end
    
    varargout{iVar} = zeros([size(varargin{iVar},1) nK]);
end

for iK = 1:1:nK
    indicesForK = k(iK)-dk/2 <= self.Kh(1,:) & self.Kh(1,:) < k(iK)+dk/2;
    for iVar=1:length(varargin)
        for iMode = 1:size(varargin{iVar},1)
            varargout{iVar}(iMode,iK) =  sum(varargin{iVar}(iMode,indicesForK));
        end
    end
end

end