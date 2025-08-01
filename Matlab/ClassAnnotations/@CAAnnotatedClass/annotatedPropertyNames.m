function names = annotatedPropertyNames(self)
% retrieve the names of all available variables
%
% - Topic: Utility function — Metadata
arguments (Input)
    self CAAnnotatedClass {mustBeNonempty}
end
arguments (Output)
    names string
end
names = self.propertyAnnotationNameMap.keys;
end