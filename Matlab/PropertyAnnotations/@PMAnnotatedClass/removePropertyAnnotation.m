function removePropertyAnnotation(self,variableAnnotation)
% add a variable annotation
%
% - Topic: Utility function — Metadata
arguments
    self PMAnnotatedClass {mustBeNonempty}
    variableAnnotation (1,:) PMPropertyAnnotation {mustBeNonempty}
end
for i=1:length(variableAnnotation)
    self.propertyAnnotationNameMap{variableAnnotation(i).name} = [];
end
end