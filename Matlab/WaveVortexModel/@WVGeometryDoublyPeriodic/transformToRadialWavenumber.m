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

% % Thi is the final output axis for wavenumber
% k = self.kRadial;
% dk = k(2)-k(1);
% 
% nK = length(k);
% Kh = self.Kh;
% 
% varargout = cell(size(varargin));
% spectralMatrixSize = self.spectralMatrixSize;
% for iVar=1:length(varargin)
%     if size(varargin{iVar},2) ~= spectralMatrixSize(2)
%         error('The input matrix must be of size [Nj Nkl]');
%     end
% 
%     varargout{iVar} = zeros([size(varargin{iVar},1) nK]);
% end
% 
% for iK = 1:1:nK
%     indicesForK = k(iK)-dk/2 <= Kh(1,:) & Kh(1,:) < k(iK)+dk/2;
%     for iVar=1:length(varargin)
%         for iMode = 1:size(varargin{iVar},1)
%             varargout{iVar}(iMode,iK) =  sum(varargin{iVar}(iMode,indicesForK));
%         end
%     end
% end

% ---- assume these already exist in your object ----
% kRadial: 1×nK vector of bin centers (monotone)
% Kh:      1×nCols vector of horizontal wavenumber magnitudes for each (k,l) column
k  = self.kRadial(:).';         % 1×nK
Kh = self.Kh(1,:);              % 1×nCols (matches your original usage)

nK    = numel(k);
nCols = numel(Kh);

% Build robust bin edges (works for nonuniform k)
if nK == 1
    % Single bin: catch everything
    dk     = max(1, 0.5*k(1));
    edges  = [k(1)-dk, k(1)+dk];
else
    mid    = 0.5*(k(1:end-1) + k(2:end));
    left   = [k(1) - (k(2)-k(1))/2, mid];
    right  = [mid, k(end) + (k(end)-k(end-1))/2];
    edges  = [left(1), mid, right(end)];
    % Compact form:
    edges  = [-Inf, mid, +Inf];         % If you prefer open-ended ends
end

% Assign each Kh to a radial bin
bin = discretize(Kh, edges);            % 1..nK or NaN for out-of-range
valid = ~isnan(bin);

% Build sparse binning matrix S: (nCols × nK), S(c, b) = 1 if column c is in bin b
S = sparse(find(valid), bin(valid), 1, nCols, nK, nnz(valid));

% Now sum columns by radial bin for each input in one shot: (nModes×nCols) * (nCols×nK)
nVars = numel(varargin);
varargout = cell(1, nVars);
for iVar = 1:nVars
    A = varargin{iVar};                 % size: nModes × nCols
    % Safety checks (optional but helpful during refactors)
    % assert(size(A,2) == nCols, 'Input %d has wrong number of columns.', iVar);

    % Matrix multiply does the grouped sum; supports real/complex
    varargout{iVar} = A * S;            % nModes × nK
end

end