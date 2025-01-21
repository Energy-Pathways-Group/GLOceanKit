function val = variableAnnotationWithName(self,name)
% retrieve a WVVariableAnnotation by name
%
% - Topic: Utility function — Metadata
arguments
    self PMAnnotatedClass {mustBeNonempty}
    name char {mustBeNonempty}
end
val = self.variableAnnotationNameMap{name};
end