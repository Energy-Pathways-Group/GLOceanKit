function val = dimensionAnnotationWithName(self,name)
% retrieve a WVDimension by name
%
% Dimension annotations are used to described the coordinate dimensions
% used in the WVTransform. These annotations are used to add variables to
% NetCDF files, and also used to generate the online documentation.
%
% Usage:
%
% ```matlab
% dimension = wvt.dimensionAnnotationWithName('x');
% ```
%
% - Topic: Metadata — Dimensions
% - Declaration: dimensionAnnotationWithName(name)
% - Parameter name: string of dimension name
% - Returns dimension: object of WVDimensionAnnotation type
arguments (Input)
    self CAAnnotatedClass {mustBeNonempty}
    name char {mustBeNonempty}
end
arguments (Output)
    val CADimensionProperty
end
val = self.dimensionAnnotationNameMap{name};
end