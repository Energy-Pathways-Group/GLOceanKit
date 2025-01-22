function addPropertyAnnotation(self,variableAnnotation)
% add a variable annotation
%
% - Topic: Utility function — Metadata
arguments
    self PMAnnotatedClass {mustBeNonempty}
    variableAnnotation (1,:) PMPropertyAnnotation
end
for i=1:length(variableAnnotation)
    self.propertyAnnotationNameMap{variableAnnotation(i).name} = variableAnnotation(i);
end
end