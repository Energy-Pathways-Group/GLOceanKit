function dimensions = defaultTransformDimensions()
dimensions = TransformDimension.empty(0,0);

dimensions(end+1) = TransformDimension('x', 'm', 'x-coordinate dimension');
dimensions(end+1) = TransformDimension('y', 'm', 'y-coordinate dimension');
dimensions(end+1) = TransformDimension('z', 'm', 'z-coordinate dimension');
dimensions(end+1) = TransformDimension('k', 'radians/m', 'wavenumber-coordinate dimension in the x-direction');
dimensions(end+1) = TransformDimension('l', 'radians/m', 'wavenumber-coordinate dimension in the y-direction');
dimensions(end+1) = TransformDimension('j', 'mode number', 'vertical mode number');
dimensions(end+1) = TransformDimension('kRadial', 'radians/m', 'isotropic wavenumber dimension');

end