function val = propertyAnnotationWithName(self,name)
% retrieve a WVPropertyAnnotation by name
%
% - Topic: Utility function — Metadata
arguments
    self WVTransform {mustBeNonempty}
    name char {mustBeNonempty}
end
val = self.propertyAnnotationNameMap(name);
end