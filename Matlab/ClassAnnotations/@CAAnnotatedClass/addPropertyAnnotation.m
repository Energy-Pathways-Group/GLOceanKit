function addPropertyAnnotation(self,propertyAnnotation)
% add a variable annotation
%
% - Topic: Utility function — Metadata
arguments
    self CAAnnotatedClass {mustBeNonempty}
    propertyAnnotation (1,:) CAPropertyAnnotation
end
for i=1:length(propertyAnnotation)
    if isa(propertyAnnotation(i),'CANumericProperty')
        isKnownDimension = isKey(self.dimensionAnnotationNameMap,propertyAnnotation(i).dimensions);
        if any(~isKnownDimension)
            error("Unable to find dimension (%s) for the numeric property %s",strjoin(string(propertyAnnotation(i).dimensions(~isKnownDimension)),","),propertyAnnotation(i).name);
        end
    end
    self.propertyAnnotationNameMap{propertyAnnotation(i).name} = propertyAnnotation(i);
    if isa(propertyAnnotation(i),'CADimensionProperty')
        self.dimensionAnnotationNameMap{propertyAnnotation(i).name} = propertyAnnotation(i);
    end
end
notify(self,'propertyAnnotationsDidChange')
end