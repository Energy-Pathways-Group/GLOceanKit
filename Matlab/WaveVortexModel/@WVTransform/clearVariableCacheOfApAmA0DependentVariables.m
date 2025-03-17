function clearVariableCacheOfApAmA0DependentVariables(self)
% clear the internal cache
%
% - Topic: Internal
%self.variableCache(intersect(self.variableCache.keys,self.wvCoefficientDependentVariablesNameMap.keys)) = [];
keys = intersect(self.variableCache.keys,self.wvCoefficientDependentVariablesNameMap.keys);
if ~isempty(keys)
    self.variableCache.remove(cellstr(keys));
end
end