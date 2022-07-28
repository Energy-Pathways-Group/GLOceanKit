function methodAnnotations = defaultMethodAnnotations()
% return array of WVAnnotations to annotate the methods
%
% This function lets us efficiently annotate all the wave vortex transform
% methods
%
% - Topic: Internal
% - Declaration: methodAnnotations = defaultMethodAnnotations()
% - Returns methodAnnotations: array of WVAnnotations instances
methodAnnotations = WVAnnotation.empty(0,0);

methodAnnotations(end+1) = WVAnnotation('transformFromSpatialDomainWithF', 'transforms from the spatial domain (x,y,z) to the spectral domain (k,l,j) using the F-modes');
methodAnnotations(end+1) = WVAnnotation('transformFromSpatialDomainWithG', 'transforms from the spatial domain (x,y,z) to the spectral domain (k,l,j) using the G-modes');
methodAnnotations(end+1) = WVAnnotation('transformToSpatialDomainWithF', 'transforms from the spectral domain (k,l,j) to the spatial domain (x,y,z) using the F-modes');
methodAnnotations(end+1) = WVAnnotation('transformToSpatialDomainWithG', 'transforms from the spectral domain (k,l,j) to the spatial domain (x,y,z) using the G-modes');
methodAnnotations(end+1) = WVAnnotation('transformToSpatialDomainWithFAllDerivatives', 'transforms from the spectral domain (k,l,j) to the spatial domain (x,y,z) using the F-modes, returning the transformed variable an its derivatives.');
methodAnnotations(end+1) = WVAnnotation('transformToSpatialDomainWithGAllDerivatives', 'transforms from the spectral domain (k,l,j) to the spatial domain (x,y,z) using the G-modes, returning the transformed variable an its derivatives.');
methodAnnotations(end+1) = WVAnnotation('diffZF', 'differentiates a variable of (x,y,z) by projecting onto the F-modes, differentiating, and transforming back to (x,y,z)');
methodAnnotations(end+1) = WVAnnotation('diffZG', 'differentiates a variable of (x,y,z) by projecting onto the G-modes, differentiating, and transforming back to (x,y,z)');
methodAnnotations(end+1) = WVAnnotation('xyzGrid', 'returns the X, Y, Z coordinate matrices',detailedDescription='- topic: Domain Attributes — Grid');
methodAnnotations(end+1) = WVAnnotation('kljGrid', 'returns the K, L, J coordinate matrices',detailedDescription='- topic: Domain Attributes — Grid');
end