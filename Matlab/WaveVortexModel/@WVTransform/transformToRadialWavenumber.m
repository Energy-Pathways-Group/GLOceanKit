function [varargout] = transformToRadialWavenumber(self,varargin)     
% transforms in the spectral domain from (k,l,j) to (kRadial,j)
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
% nexttile, pcolor(wvt.j,wvt.kRadial,Ekj), shading flat
% nexttile, plot(wvt.j,sum(Ekj,1))
% ```
%
% - Topic: Operations — Transformations
% - Declaration: [varargout] = transformToRadialWavenumber(varargin) 
% - Parameter varargin: variables with dimensions $$(k,l)$$ or $$(k,l,j)$$
% - Returns varargout: variables with dimensions $$(kRadial)$$ or $$(kRadial,j)$$
Kh = self.Kh;

% Thi is the final output axis for wavenumber
k = self.kRadial;
dk = k(2)-k(1);

RedundantCoefficients = WVTransform.redundantHermitianCoefficients(Kh);
OmNyquist = self.maskForNyquistModes();
nK = length(k);

varargout = cell(size(varargin));
for iVar=1:length(varargin)
    thesize = size(varargin{iVar});
    if length(thesize) == 2
        newsize = [nK 1];
    else
        newsize = cat(2,nK,thesize(3:end));
    end
    varargout{iVar} = zeros(newsize);
end

for iK = 1:1:nK
    indicesForK = find( k(iK)-dk/2 <= squeeze(Kh(:,:,1)) & squeeze(Kh(:,:,1)) < k(iK)+dk/2  & ~squeeze(OmNyquist(:,:,1)) & ~squeeze(RedundantCoefficients(:,:,1))  );
    for iIndex = 1:length(indicesForK)
        [i,m] = ind2sub([size(Kh,1) size(Kh,2)],indicesForK(iIndex));
        if i+m==2
            prefactor = 1;
        else
            prefactor = 2;
        end
        for iMode = 1:length(self.j)
            for iVar=1:length(varargin)
                if ismatrix(varargin{iVar})
                    continue;
                end
                varargout{iVar}(iK,iMode) = varargout{iVar}(iK,iMode) + prefactor*varargin{iVar}(i,m,iMode);
            end
        end
        for iVar=1:length(varargin)
            if ~ismatrix(varargin{iVar})
                continue;
            end
            varargout{iVar}(iK) = varargout{iVar}(iK) + prefactor*varargin{iVar}(i,m);
        end
    end
end
end