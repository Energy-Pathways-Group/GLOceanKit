function names = variableNames(self)
% retrieve the names of all available variables
%
% - Topic: Utility function — Metadata
arguments
    self PMAnnotatedClass {mustBeNonempty}
end
names = self.variableAnnotationNameMap.keys;
end