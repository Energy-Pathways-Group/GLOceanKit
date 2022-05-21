function [varargout] = VariablesAtTimePosition(self,t,x,y,z,interpolationMethod,varargin)
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
    varargout = cell(size(varargin));
    [varargout{:}] = self.InternalVariablesAtTimePosition(t,x,y,z,interpolationMethod,varargin{:});
    if ~isempty(self.offgridModes.k_ext)
        varargoutExt = cell(size(varargin));
        [varargoutExt{:}] = self.ExternalVariablesAtTimePosition(t,x,y,z,varargin{:});
        for iArg=1:length(varargout)
            varargout{iArg} = varargout{iArg} + varargoutExt{iArg};
        end
    end
end
