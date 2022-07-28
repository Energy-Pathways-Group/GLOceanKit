function dimensions = defaultDimensionAnnotations()
% return array of TransformDimensions initialized by default
%
% This function lets us efficiently annotate all the coordinate dimensions.
%
% - Topic: Internal
% - Declaration: dimensions = defaultDimensionAnnotations()
% - Returns dimensions: array of WVDimensionAnnotation instances
dimensions = WVDimensionAnnotation.empty(0,0);

dimensions(end+1) = WVDimensionAnnotation('t', 's', 'time dimension');
dimensions(end+1) = WVDimensionAnnotation('x', 'm', 'x-coordinate dimension');
dimensions(end+1) = WVDimensionAnnotation('y', 'm', 'y-coordinate dimension');
dimensions(end+1) = WVDimensionAnnotation('z', 'm', 'z-coordinate dimension');
dimensions(end+1) = WVDimensionAnnotation('k', 'rad/m', 'wavenumber-coordinate dimension in the x-direction');
dimensions(end+1) = WVDimensionAnnotation('l', 'rad/m', 'wavenumber-coordinate dimension in the y-direction');
dimensions(end+1) = WVDimensionAnnotation('j', 'mode number', 'vertical mode number');
dimensions(end+1) = WVDimensionAnnotation('kRadial', 'rad/m', 'isotropic wavenumber dimension');

end