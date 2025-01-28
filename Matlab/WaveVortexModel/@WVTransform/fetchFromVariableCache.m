function varargout = fetchFromVariableCache(self,varargin)
% retrieve a set of variables from the internal cache
%
% - Topic: Internal
varargout = cell(size(varargin));
for iVar=1:length(varargin)
    if isKey(self.variableCache,varargin{iVar})
        varargout{iVar} = self.variableCache{varargin{iVar}};
    else
        varargout{iVar} = [];
    end
end
end