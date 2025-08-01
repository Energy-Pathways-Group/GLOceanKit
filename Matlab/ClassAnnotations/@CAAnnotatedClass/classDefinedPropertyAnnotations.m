function propertyAnnotations = classDefinedPropertyAnnotations()
% return array of WVPropertyAnnotation initialized by default
%
% This function returns annotations for all properties of the
% CAAnnotatedClass class.
%
% - Topic: Developer
% - Declaration: propertyAnnotations = CAAnnotatedClass.classVariableAnnotations()
% - Returns propertyAnnotations: array of WVPropertyAnnotation instances
arguments (Output)
    propertyAnnotations CAPropertyAnnotation
end
propertyAnnotations = CAPropertyAnnotation.empty(0,0);
end