function propertyAnnotations = defaultPropertyAnnotations()
% return array of WVPropertyAnnotation initialized by default
%
% This function lets us efficiently annotate all the wave vortex transform
% properties
%
% - Topic: Internal
% - Declaration: propertyAnnotations = defaultPropertyAnnotations()
% - Returns propertyAnnotations: array of WVPropertyAnnotation instances
propertyAnnotations = WVPropertyAnnotation.empty(0,0);

propertyAnnotations(end+1) = WVPropertyAnnotation('F',{'k','l','j'},'', 'normalization factor for the $$F$$ modes',isComplex=0);
propertyAnnotations(end+1) = WVPropertyAnnotation('G',{'k','l','j'},'', 'normalization factor for the $$G$$ modes.',isComplex=0);
propertyAnnotations(end+1) = WVPropertyAnnotation('h',{'k','l','j'},'m', 'equivalent depth of each mode');

propertyAnnotations(end+1) = WVPropertyAnnotation('cg_x',{'k','l','j'},'m/s', 'wave group speed in the x-direction', detailedDescription='- topic: Domain Attributes');
propertyAnnotations(end+1) = WVPropertyAnnotation('cg_y',{'k','l','j'},'m/s', 'wave group speed in the y-direction', detailedDescription='- topic: Domain Attributes');
propertyAnnotations(end+1) = WVPropertyAnnotation('cg_z',{'k','l','j'},'m/s', 'wave group speed in the z-direction', detailedDescription='- topic: Domain Attributes');



end