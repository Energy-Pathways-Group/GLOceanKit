function clearVariableCache(self)
% clear the internal cache
%
% - Topic: Internal
self.variableCache = configureDictionary("string","cell");
end