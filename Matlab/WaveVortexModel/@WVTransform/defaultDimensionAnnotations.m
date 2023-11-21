function dimensions = defaultDimensionAnnotations()
% return array of TransformDimensions initialized by default
%
% This function lets us efficiently annotate all the coordinate dimensions.
%
% - Topic: Internal
% - Declaration: dimensions = defaultDimensionAnnotations()
% - Returns dimensions: array of WVDimensionAnnotation instances
dimensions = WVDimensionAnnotation.empty(0,0);

dimensions(end+1) = WVDimensionAnnotation('t', 's', 'time coordinate');
dimensions(end).attributes('standard_name') = 'time';
dimensions(end).attributes('axis') = 'T';

dimensions(end+1) = WVDimensionAnnotation('x', 'm', 'x coordinate');
dimensions(end).attributes('standard_name') = 'projection_x_coordinate';
dimensions(end).attributes('axis') = 'X';

dimensions(end+1) = WVDimensionAnnotation('y', 'm', 'y coordinate');
dimensions(end).attributes('standard_name') = 'projection_y_coordinate';
dimensions(end).attributes('axis') = 'Y';

dimensions(end+1) = WVDimensionAnnotation('z', 'm', 'z coordinate');
dimensions(end).attributes('standard_name') = 'height_above_mean_sea_level';
dimensions(end).attributes('positive') = 'up';
dimensions(end).attributes('axis') = 'Z';

dimensions(end+1) = WVDimensionAnnotation('k', 'rad/m', 'wavenumber coordinate in the x-direction');
dimensions(end+1) = WVDimensionAnnotation('l', 'rad/m', 'wavenumber coordinate in the y-direction');
dimensions(end+1) = WVDimensionAnnotation('j', 'mode number', 'vertical mode number');
dimensions(end+1) = WVDimensionAnnotation('kRadial', 'rad/m', 'isotropic wavenumber dimension');

end