function propertyAnnotations = propertyAnnotationForKnownVariable(variableName,options)
% This is one of two functions that returns operations for computing
% standard variables.
arguments (Input,Repeating)
    variableName char
end
arguments (Input)
    options.spectralDimensionNames = {'j','kl'}
    options.spatialDimensionNames = {'x','y','z'}
    options.primaryFlowComponents WVFlowComponent = WVFlowComponent.empty(0,0)
    options.flowComponent WVFlowComponent = WVFlowComponent.empty(0,0)
end
arguments (Output)
    propertyAnnotations CAPropertyAnnotation
end

operations = WVTransform.classDefinedOperationForKnownVariable(variableName{:},spectralDimensionNames=options.spectralDimensionNames,spatialDimensionNames=options.spatialDimensionNames);
propertyAnnotations = CAPropertyAnnotation.empty(0,0);
for iOp = 1:length(operations)
    propertyAnnotations(iOp) = operations(iOp).outputVariables;
end
end