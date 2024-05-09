function removeVariableAnnotations(self,variableAnnotation)
% add a variable annotation
%
% - Topic: Utility function — Metadata
arguments
    self WVTransform {mustBeNonempty}
    variableAnnotation (1,:) WVVariableAnnotation {mustBeNonempty}
end
for i=1:length(variableAnnotation)
    remove(self.variableAnnotationNameMap,variableAnnotation(i).name);
    if isKey(self.timeDependentVariablesNameMap,variableAnnotation(i).name)
        remove(self.timeDependentVariablesNameMap,variableAnnotation(i).name);
    end
end
end