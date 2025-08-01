function [varargout] = variableWithName(self, variableNames)
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
didFetchAll = 0;
while didFetchAll ~= 1
    [varargout{:}] = self.fetchFromVariableCache(variableNames{:});

    for iVar=1:length(varargout)
        if isempty(varargout{iVar})
            % now go compute it, and then try fetching from
            % cach
            modelVar = self.operationVariableNameMap(variableNames{iVar});
            self.performOperationWithName(modelVar.modelOp.name);
            didFetchAll = 0;
            break;
        else
            didFetchAll = 1;
        end
    end
end

% % External variables may not all be supported.
% if ~isempty(self.offgridModes) && ~isempty(self.offgridModes.k_ext)
%     varargoutExt = cell(size(variableNames));
%     [varargoutExt{:}] = self.externalVariableFieldsAtTime(self.t, variableNames{:});
%     for iArg=1:length(varargout)
%         varargout{iArg} = varargout{iArg} + varargoutExt{iArg};
%     end
% end
end