function [varargout] = transformToKLAxes(self,varargin)     
% transforms in the spectral domain from (j,kl) to (kAxis,lAxis,j)
%
%
% The following example takes the total energy of the wave part of
% flow, transforms it to the (kAxis,lAxis,j) grid, then sums all the energy
% along the mode (j) dimension.
%
% ```matlab
% figure
% TE = wvt.Apm_TE_factor .* (abs(wvt.Ap).^2 + abs(wvt.Am).^2);
% pcolor(wvt.kAxis,wvt.lAxis,log10(sum(wvt.transformToKLAxes(TE),3)).'), shading flat
% ```
%
% - Topic: Operations — Transformations
% - Declaration: [varargout] = transformToKLAxes(varargin) 
% - Parameter varargin: variables with dimensions $$(j,kl)$$
% - Returns varargout: variables with dimensions (kAxis,lAxis,j)

varargout = cell(size(varargin));
for iVar=1:length(varargin)
    if size(varargin{iVar},2) ~= self.Nkl && numel(varargin{iVar}) ~= self.Nkl
        error('The input matrix must have the second dimension of length Nkl, or be a vector of that length.');
    end
    
    if size(varargin{iVar},2) == self.Nkl
    varargout{iVar} = zeros([size(varargin{iVar},1) length(self.kAxis) length(self.lAxis)]);
    else
        varargout{iVar} = [];
    end
end

for iVar=1:length(varargin)
    varargout{iVar} =  fftshift(fftshift(self.transformFromWVGridToDFTGrid(varargin{iVar}),1),2);
    % If you don't want it to swap the axis ordering:
    % for iMode=1:size(varargin{iVar},1)
    %     varargout{iVar}(iMode,:,:) =  fftshift(fftshift(self.horizontalModes.transformFromWVGridToDFTGrid(varargin{iVar}(iMode,:)),2),3);
    % end
end

end