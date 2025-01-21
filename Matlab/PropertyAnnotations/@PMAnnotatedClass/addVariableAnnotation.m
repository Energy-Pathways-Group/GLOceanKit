function addVariableAnnotation(self,variableAnnotation)
% add a variable annotation
%
% - Topic: Utility function — Metadata
arguments
    self PMAnnotatedClass {mustBeNonempty}
    variableAnnotation (1,:) PMVariableAnnotation {mustBeNonempty}
end
for i=1:length(variableAnnotation)
    self.variableAnnotationNameMap{variableAnnotation(i).name} = variableAnnotation(i);
end
end