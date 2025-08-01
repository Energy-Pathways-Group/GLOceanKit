function addToVariableCache(self,name,var)
% add variable to internal cache, in case it is needed again
%
% - Topic: Internal
self.variableCache{name} = var;
end