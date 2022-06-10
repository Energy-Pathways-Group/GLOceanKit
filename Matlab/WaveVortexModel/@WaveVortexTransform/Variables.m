function [varargout] = Variables(self, varargin)
% Primary method for accessing the dynamical variables on the
% internal grid.

varargout = cell(size(varargin));
[varargout{:}] = self.stateVariables(varargin{:});

% External variables may not all be supported.
if ~isempty(self.offgridModes) && ~isempty(self.offgridModes.k_ext)
    varargoutExt = cell(size(varargin));
    [varargoutExt{:}] = self.ExternalVariableFieldsAtTime(self.t, varargin{:});
    for iArg=1:length(varargout)
        varargout{iArg} = varargout{iArg} + varargoutExt{iArg};
    end
end
end