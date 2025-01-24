function addPropertyAnnotation(self,propertyAnnotation)
% add a variable annotation
%
% - Topic: Utility function — Metadata
arguments
    self CAAnnotatedClass {mustBeNonempty}
    propertyAnnotation (1,:) CAPropertyAnnotation
end
for i=1:length(propertyAnnotation)
    self.propertyAnnotationNameMap{propertyAnnotation(i).name} = propertyAnnotation(i);
    if isa(propertyAnnotation(i),'CADimensionProperty')
        self.dimensionAnnotationNameMap{propertyAnnotation(i).name} = propertyAnnotation(i);
    end
end
end