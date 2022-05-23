function [varargout] = Variables(self, varargin)
% Primary method for accessing the dynamical variables on the
% internal grid.
%
% Valid variable options are 'u', 'v', 'w', 'rho_prime', and
% 'zeta'.
varargout = cell(size(varargin));
[varargout{:}] = self.TransformVariables(varargin{:});
if ~isempty(self.offgridModes.k_ext)
    varargoutExt = cell(size(varargin));
    [varargoutExt{:}] = self.ExternalVariableFieldsAtTime(self.t, varargin{:});
    for iArg=1:length(varargout)
        varargout{iArg} = varargout{iArg} + varargoutExt{iArg};
    end
end
end