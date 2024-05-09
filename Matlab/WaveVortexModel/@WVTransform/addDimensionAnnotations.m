function addDimensionAnnotations(self,dimensionAnnotation)
% add one or more WVDimensions
%
% - Topic: Utility function — Metadata
arguments
    self WVTransform {mustBeNonempty}
    dimensionAnnotation (1,:) WVDimensionAnnotation {mustBeNonempty}
end
for i=1:length(dimensionAnnotation)
    self.dimensionAnnotationNameMap(dimensionAnnotation(i).name) = dimensionAnnotation(i);
end
end