function [varargout] = variables(self, variableNames)
% Primary method for accessing the dynamical variables
%
% Computes (or retrieves from cache) any known state variables.
%
% - Topic: State Variables
% - Declaration: [varargout] = variables(self, variables)
% - Parameter variables: strings of variable names.
arguments
    self WVTransform {mustBeNonempty}
end
arguments (Repeating)
    variableNames char
end

varargout = cell(size(variableNames));
[varargout{:}] = self.stateVariables(variableNames{:});

% External variables may not all be supported.
if ~isempty(self.offgridModes) && ~isempty(self.offgridModes.k_ext)
    varargoutExt = cell(size(variableNames));
    [varargoutExt{:}] = self.externalVariableFieldsAtTime(self.t, variableNames{:});
    for iArg=1:length(varargout)
        varargout{iArg} = varargout{iArg} + varargoutExt{iArg};
    end
end
end