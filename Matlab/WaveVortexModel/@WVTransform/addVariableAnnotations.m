function addVariableAnnotations(self,variableAnnotation)
% add a variable annotation
%
% - Topic: Utility function — Metadata
arguments
    self WVTransform {mustBeNonempty}
    variableAnnotation (1,:) WVVariableAnnotation {mustBeNonempty}
end
for i=1:length(variableAnnotation)
    self.variableAnnotationNameMap(variableAnnotation(i).name) = variableAnnotation(i);
    if variableAnnotation(i).isVariableWithLinearTimeStep == 1 && ~isKey(self.timeDependentVariablesNameMap,variableAnnotation(i).name)
        self.timeDependentVariablesNameMap(variableAnnotation(i).name) = variableAnnotation(i);
    end
end
end