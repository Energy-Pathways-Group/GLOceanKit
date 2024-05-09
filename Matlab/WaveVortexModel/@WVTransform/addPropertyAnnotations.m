function addPropertyAnnotations(self,propertyAnnotation)
% add a property annotation
%
% - Topic: Utility function — Metadata
arguments
    self WVTransform {mustBeNonempty}
    propertyAnnotation (1,:) WVPropertyAnnotation {mustBeNonempty}
end
for i=1:length(propertyAnnotation)
    self.propertyAnnotationNameMap(propertyAnnotation(i).name) = propertyAnnotation(i);
end
end