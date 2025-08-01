function removePropertyAnnotation(self,propertyAnnotation)
% add a variable annotation
%
% - Topic: Utility function — Metadata
arguments
    self CAAnnotatedClass {mustBeNonempty}
    propertyAnnotation (1,:) CAPropertyAnnotation {mustBeNonempty}
end
for i=1:length(propertyAnnotation)
    self.propertyAnnotationNameMap{propertyAnnotation(i).name} = [];
    if isa(propertyAnnotation(i),'CADimensionProperty')
        self.dimensionAnnotationNameMap{propertyAnnotation(i).name} = [];
    end
end
notify(self,'propertyAnnotationsDidChange')
end