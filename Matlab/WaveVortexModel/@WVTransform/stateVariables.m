function [varargout] = stateVariables(self, varargin)
% retrieve variables either from cache or by computation
%
% - Topic: Internal
varargout = cell(size(varargin));

didFetchAll = 0;
while didFetchAll ~= 1
    [varargout{:}] = self.fetchFromVariableCache(varargin{:});

    for iVar=1:length(varargout)
        if isempty(varargout{iVar})
            % now go compute it, and then try fetching from
            % cach
            modelVar = self.variableAnnotationNameMap(varargin{iVar});
            self.performOperationWithName(modelVar.modelOp.name);
            didFetchAll = 0;
            break;
        else
            didFetchAll = 1;
        end
    end
end
end