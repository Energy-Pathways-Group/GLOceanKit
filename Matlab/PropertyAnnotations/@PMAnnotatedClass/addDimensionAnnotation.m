function addDimensionAnnotation(self,dimensionAnnotation)
% add one or more PMDimensionAnnotation
%
% Dimension annotations are used to described the coordinate dimensions
% used in the WVTransform. These annotations are used to add variables to
% NetCDF files, and also used to generate the online documentation.
%
% In general, users should not need to add dimensions. The default
% dimensions are added by the WVTransform upon initialization.
% Usage:
%
% ```matlab
% dimensions = WVDimensionAnnotation('x', 'm', 'x coordinate');
% dimensions.attributes('standard_name') = 'projection_x_coordinate';
% dimensions.attributes('axis') = 'X';
%
% wvt.addDimensionAnnotations(dimension);
% ```
% - Topic: Metadata — Dimensions
% - Declaration: addDimensionAnnotations(dimensionAnnotation)
% - Parameter dimensionAnnotation: one or more WVDimensionAnnotation objects
arguments
    self PMAnnotatedClass {mustBeNonempty}
    dimensionAnnotation (1,:) PMDimensionAnnotation {mustBeNonempty}
end
for i=1:length(dimensionAnnotation)
    self.dimensionAnnotationNameMap{dimensionAnnotation(i).name} = dimensionAnnotation(i);
end
end