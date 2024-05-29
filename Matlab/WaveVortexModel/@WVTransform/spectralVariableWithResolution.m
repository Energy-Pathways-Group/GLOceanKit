function [varargout] = spectralVariableWithResolution(self,wvtX2,varargin)
% create a new variable with different resolution
%
% Given a variable with dimensions [Nj Nkl], this returns a new variable
% with dimensions matching wvtX2. This can be either increased or decreased
% resolution.
%
% - Topic: Utility function
% - Declaration: varX2 = spectralVariableWithResolution(var,Nklj)
% - Parameter var: a variable with dimensions [Nj Nkl]
% - Parameter wvtX2: a WVTransform of different size.
% - Returns varX2: matrix the size Nklj

if ~isequal(self.Lx,wvtX2.Lx) || ~isequal(self.Ly,wvtX2.Ly) || ~isequal(self.Lz,wvtX2.Lz)
    error('These transforms are not compatible.')
end

Nkl = min(self.Nkl,wvtX2.Nkl);
Nj = min(self.Nj,wvtX2.Nj);

varargout = cell(size(varargin));
for iVar=1:length(varargin)
    if isempty(varargin{iVar})
        varargout{iVar} = [];
    else
        varargout{iVar} = zeros(wvtX2.spectralMatrixSize);
        varargout{iVar}(1:Nj,1:Nkl) = varargin{iVar}(1:Nj,1:Nkl);
    end
end

end