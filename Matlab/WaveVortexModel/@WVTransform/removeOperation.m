function removeOperation(self,transformOperation)
% remove an existing WVOperation
%
% - Topic: Utility function — Metadata
arguments
    self WVTransform {mustBeNonempty}
    transformOperation (1,1) WVOperation {mustBeNonempty}
end
self.removeVariableAnnotations(transformOperation.outputVariables);
for iVar=1:transformOperation.nVarOut
    self.operationVariableNameMap(transformOperation.outputVariables(iVar).name) = [];
end
self.operationNameMap(transformOperation.name) = [];
self.clearVariableCache();
end