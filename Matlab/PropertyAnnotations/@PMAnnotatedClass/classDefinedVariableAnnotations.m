function variableAnnotations = classDefinedVariableAnnotations()
% return array of WVPropertyAnnotation initialized by default
%
% This function returns annotations for all properties of the
% PMAnnotatedClass class.
%
% - Topic: Developer
% - Declaration: propertyAnnotations = PMAnnotatedClass.classVariableAnnotations()
% - Returns propertyAnnotations: array of WVPropertyAnnotation instances
variableAnnotations = PMVariableAnnotation.empty(0,0);
end