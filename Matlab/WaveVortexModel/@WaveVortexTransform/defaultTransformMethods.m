function transformMethods = defaultTransformMethods()
% return array of TransformAnnotations to annotate the methods
%
% This function lets us efficiently annotate all the wave vortex transform
% methods
%
% - Topic: Internal
% - Declaration: transformMethods = defaultTransformMethods()
% - Returns transformProperties: array of TransformAnnotations instances
transformMethods = TransformAnnotation.empty(0,0);

transformMethods(end+1) = TransformAnnotation('transformFromSpatialDomainWithF', 'transforms from the spatial domain (x,y,z) to the spectral domain (k,l,j) using the F-modes');
transformMethods(end+1) = TransformAnnotation('transformFromSpatialDomainWithG', 'transforms from the spatial domain (x,y,z) to the spectral domain (k,l,j) using the G-modes');
transformMethods(end+1) = TransformAnnotation('transformToSpatialDomainWithF', 'transforms from the spectral domain (k,l,j) to the spatial domain (x,y,z) using the F-modes');
transformMethods(end+1) = TransformAnnotation('transformToSpatialDomainWithG', 'transforms from the spectral domain (k,l,j) to the spatial domain (x,y,z) using the G-modes');
transformMethods(end+1) = TransformAnnotation('transformToSpatialDomainWithFAllDerivatives', 'transforms from the spectral domain (k,l,j) to the spatial domain (x,y,z) using the F-modes, returning the transformed variable an its derivatives.');
transformMethods(end+1) = TransformAnnotation('transformToSpatialDomainWithGAllDerivatives', 'transforms from the spectral domain (k,l,j) to the spatial domain (x,y,z) using the G-modes, returning the transformed variable an its derivatives.');
transformMethods(end+1) = TransformAnnotation('diffZF', 'differentiates a variable of (x,y,z) by projecting onto the F-modes, differentiating, and transforming back to (x,y,z)');
transformMethods(end+1) = TransformAnnotation('diffZG', 'differentiates a variable of (x,y,z) by projecting onto the G-modes, differentiating, and transforming back to (x,y,z)');
end