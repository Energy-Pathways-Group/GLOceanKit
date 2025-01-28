function [varargout] = variableAtPositionWithName(self,x,y,z,variableNames,options)
% Primary method for accessing the dynamical variables on the at any
% position in the domain.
%
% Computes (or retrieves from cache) any known state variables and computes
% their values at the requested positions (x,y,z)
%
% The method argument specifies how off-grid values should be interpolated.
% Use 'exact' for the slow, but accurate, spectral interpolation. Otherwise
% use 'spline' or some other method used by Matlab's interp function.
%
% - Topic: State Variables
% - Declaration: [varargout] = variableAtPositionWithName(self,x,y,z,variableNames,options)
% - Parameter x: array of x-positions
% - Parameter y: array of y-positions
% - Parameter z: array of z-positions
% - Parameter variableNames: strings of variable names.
% - Parameter interpolationMethod: (optional) `linear`,`spline`,`exact`. Default `linear`.
arguments
    self WVTransform {mustBeNonempty}
    x (1,:) double
    y (1,:) double
    z (1,:) double
end
arguments (Repeating)
    variableNames char
end
arguments
    options.interpolationMethod char {mustBeMember(options.interpolationMethod,["linear","spline","exact","finufft"])} = "linear"
end

varargout = cell(size(variableNames));
if strcmp(options.interpolationMethod,"exact")
    error('I ripped this out.')
    if isempty(self.ongridModes)
        self.ongridModes = WVOffGridTransform(self.offgridModes.verticalModes,self.offgridModes.latitude,self.offgridModes.N2Function);
        [omega, alpha, ~, ~, mode, phi, A, norm] = self.waveModesFromWaveCoefficients();
        self.ongridModes.setExternalWavesWithFrequencies(omega, alpha, mode, phi, A, norm);
    end
    [varargout{:}] = self.ongridModes.externalVariablesAtTimePosition(t,x,y,z, variableNames{:});
else
    [varargout{:}] = self.variableWithName(variableNames{:});
    [varargout{:}] = self.interpolatedFieldAtPosition(x,y,z,options.interpolationMethod,varargout{:});
end

% if ~isempty(self.offgridModes) && ~isempty(self.offgridModes.k_ext)
%     varargoutExt = cell(size(variableNames));
%     [varargoutExt{:}] = self.externalVariablesAtTimePosition(self.t,x,y,z,variableNames{:});
%     for iArg=1:length(varargout)
%         varargout{iArg} = varargout{iArg} + varargoutExt{iArg};
%     end
% end
end
