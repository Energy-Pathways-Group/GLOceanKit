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

methodAnnotations(end+1) = WVAnnotation('transformFromSpatialDomainWithFio', 'transforms from the spatial domain (z,:,:) to the spectral domain (j,:,:) using the inertial oscillation F-modes');
methodAnnotations(end+1) = WVAnnotation('transformFromSpatialDomainWithFg', 'transforms from the spatial domain (z,:,:) to the spectral domain (j,:,:) using the geostrophic F-modes');
methodAnnotations(end+1) = WVAnnotation('transformFromSpatialDomainWithGg', 'transforms from the spatial domain (z,:,:) to the spectral domain (j,:,:) using the geostrophic G-modes');
methodAnnotations(end+1) = WVAnnotation('transformToSpatialDomainWithF', 'transforms from the spectral domain (k,l,j) to the spatial domain (x,y,z) using the F-modes');
methodAnnotations(end+1) = WVAnnotation('transformToSpatialDomainWithG', 'transforms from the spectral domain (k,l,j) to the spatial domain (x,y,z) using the G-modes');
methodAnnotations(end+1) = WVAnnotation('transformToSpatialDomainWithFAllDerivatives', 'transforms from the spectral domain (k,l,j) to the spatial domain (x,y,z) using the F-modes, returning the transformed variable an its derivatives.');
methodAnnotations(end+1) = WVAnnotation('transformToSpatialDomainWithGAllDerivatives', 'transforms from the spectral domain (k,l,j) to the spatial domain (x,y,z) using the G-modes, returning the transformed variable an its derivatives.');
methodAnnotations(end+1) = WVAnnotation('diffZF', 'differentiates a variable of (x,y,z) by projecting onto the F-modes, differentiating, and transforming back to (x,y,z)');
methodAnnotations(end+1) = WVAnnotation('diffZG', 'differentiates a variable of (x,y,z) by projecting onto the G-modes, differentiating, and transforming back to (x,y,z)');
methodAnnotations(end+1) = WVAnnotation('xyzGrid', 'returns the X, Y, Z coordinate matrices',detailedDescription='- topic: Domain Attributes — Grid — Spatial');
methodAnnotations(end+1) = WVAnnotation('kljGrid', 'returns the K, L, J coordinate matrices',detailedDescription='- topic: Domain Attributes — Grid — Spectral');
methodAnnotations(end+1) = WVAnnotation('spatialMatrixSize', 'returns the size of all real-valued field variables',detailedDescription='- topic: Domain Attributes — Grid — Spatial');
methodAnnotations(end+1) = WVAnnotation('spectralMatrixSize', 'returns the size of any spectral matrix, e.g., Ap, Am, A0',detailedDescription='- topic: Domain Attributes — Grid — Spectral');
methodAnnotations(end+1) = WVAnnotation('variables', 'access the dynamical variables');
methodAnnotations(end+1) = WVAnnotation('variableAtPositionWithName', 'access the dynamical variables at any position in the domain');
end