function [varargout] = variablesAtPosition(self,x,y,z,variableNames,options)
arguments
    self WaveVortexTransform {mustBeNonempty}
    x (1,:) double
    y (1,:) double
    z (1,:) double
end
arguments (Repeating)
    variableNames char
end
arguments
    options.InterpolationMethod char {mustBeMember(options.InterpolationMethod,["linear","spline","exact"])} = "linear"
end
    % Primary method for accessing the dynamical variables on the
    % at any position or time.
    %
    % The method argument specifies how off-grid values should be
    % interpolated. Use 'exact' for the slow, but accurate,
    % spectral interpolation. Otherwise use 'spline' or some other
    % method used by Matlab's interp function.
    %
    % Valid variable options are 'u', 'v', 'w', 'rho_prime', and
    % 'zeta'.
    varargout = cell(size(variableNames));
    if strcmp(options.InterpolationMethod,"exact")
        if isempty(self.ongridModes)
            self.ongridModes = WaveVortexModelOffGrid(self.offgridModes.internalModes,self.offgridModes.latitude,self.offgridModes.N2Function);
            [omega, alpha, ~, ~, mode, phi, A, norm] = self.waveModesFromWaveCoefficients();
            self.ongridModes.setExternalWavesWithFrequencies(omega, alpha, mode, phi, A, norm);
        end
        [varargout{:}] = self.ongridModes.externalVariablesAtTimePosition(t,x,y,z, variableNames{:}); 
    else
        [varargout{:}] = self.stateVariables(variableNames{:});
        [varargout{:}] = self.interpolatedFieldAtPosition(x,y,z,options.InterpolationMethod,varargout{:});
    end

    if ~isempty(self.offgridModes) && ~isempty(self.offgridModes.k_ext)
        varargoutExt = cell(size(variableNames));
        [varargoutExt{:}] = self.externalVariablesAtTimePosition(self.t,x,y,z,variableNames{:});
        for iArg=1:length(varargout)
            varargout{iArg} = varargout{iArg} + varargoutExt{iArg};
        end
    end
end
