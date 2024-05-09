function clearVariableCacheOfTimeDependentVariables(self)
% clear the internal cache of variables that claim to be time dependent
%
% - Topic: Internal
remove(self.variableCache,intersect(self.variableCache.keys,self.timeDependentVariablesNameMap.keys));
end