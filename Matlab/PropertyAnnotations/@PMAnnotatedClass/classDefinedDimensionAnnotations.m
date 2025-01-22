function dimensions = classDefinedDimensionAnnotations()
% return array of PMDimensionAnnotation to annotate the
% dimensions
%
% This function returns annotations for all dimensions of the
% PMAnnotatedClass class.
%
% - Topic: Developer
% - Declaration: dimensionAnnotations = PMAnnotatedClass.classDimensionAnnotations()
% - Returns dimensionAnnotations: array of PMDimensionAnnotation instances
arguments (Output)
    dimensions PMDimensionAnnotation
end
dimensions = PMDimensionAnnotation.empty(0,0);
end