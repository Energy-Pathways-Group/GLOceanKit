function names = propertyNames(self)
% retrieve the names of all available variables
%
% - Topic: Utility function — Metadata
arguments
    self PMAnnotatedClass {mustBeNonempty}
end
names = self.propertyAnnotationNameMap.keys;
end