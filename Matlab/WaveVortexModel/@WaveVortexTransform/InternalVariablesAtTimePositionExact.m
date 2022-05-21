%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% These functions return the exact/spectral velocity field for the
% internal gridded field
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [varargout] = InternalVariablesAtTimePositionExact(self,t,x,y,z,varargin)
% Return the velocity field associated with the gridded
% velocity field, but using spectral interpolation, rather than
% the FFT grid.
% size(x) = [N 1]
% size(phi) = [1 M]
if isempty(self.ongridModes)
    self.ongridModes = WaveVortexModelOffGrid(self.offgridModes.internalModes,self.offgridModes.latitude,self.offgridModes.N2Function);
    [omega, alpha, ~, ~, mode, phi, A, norm] = self.WaveCoefficientsFromGriddedWaves();
    self.ongridModes.SetExternalWavesWithFrequencies(omega, alpha, mode, phi, A, norm);
end

% if ~isempty(self.ongridModes.k_ext)
varargout = cell(size(varargin));
    [varargout{:}] = self.ongridModes.ExternalVariablesAtTimePosition(t,x,y,z, varargin{:});
% end
end