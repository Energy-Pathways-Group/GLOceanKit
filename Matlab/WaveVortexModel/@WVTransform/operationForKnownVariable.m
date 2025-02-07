function operations = operationForKnownVariable(self,variableName,options)
% This is one of two functions that returns operations for computing
% standard variables.
arguments (Input)
    self WVTransform
end
arguments (Input,Repeating)
    variableName char
end
arguments (Input)
    options.flowComponent WVFlowComponent = WVFlowComponent.empty(0,0)
end
arguments (Output)
    operations WVOperation
end
operations = WVTransform.classDefinedOperationForKnownVariable(variableName{:},flowComponent=options.flowComponent,spatialDimensionNames=self.spatialDimensionNames,spectralDimensionNames=self.spectralDimensionNames,totalFlowComponent=self.totalFlowComponent);
end