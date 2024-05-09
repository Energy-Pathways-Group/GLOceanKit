function names = variableNames(self)
% retrieve the names of all available variables
%
% - Topic: Utility function — Metadata
arguments
    self WVTransform {mustBeNonempty}
end
names = self.variableAnnotationNameMap.keys;
end