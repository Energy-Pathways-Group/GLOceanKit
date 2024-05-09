function val = dimensionAnnotationWithName(self,name)
% retrieve a WVDimension by name
%
% - Topic: Utility function — Metadata
arguments
    self WVTransform {mustBeNonempty}
    name char {mustBeNonempty}
end
val = self.dimensionAnnotationNameMap(name);
end