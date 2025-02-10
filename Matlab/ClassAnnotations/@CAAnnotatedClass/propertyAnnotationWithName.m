function val = propertyAnnotationWithName(self,name)
% retrieve a WVVariableAnnotation by name
%
% - Topic: Utility function — Metadata
arguments (Input)
    self CAAnnotatedClass {mustBeNonempty}
    name string
end
arguments (Output)
    val CAPropertyAnnotation
end
% Cannot call val = self.propertyAnnotationNameMap{name}; because of a bug?
val = CAPropertyAnnotation.empty(0,0);
for i=1:length(name)
    val(i) = self.propertyAnnotationNameMap{name(i)};
end

end