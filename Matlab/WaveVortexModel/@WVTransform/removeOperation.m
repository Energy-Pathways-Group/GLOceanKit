function removeOperation(self,operation)
% remove an existing WVOperation
%
% - Topic: Utility function — Metadata
arguments
    self WVTransform {mustBeNonempty}
    operation (1,1) WVOperation {mustBeNonempty}
end
self.removePropertyAnnotation(operation.outputVariables);
for iVar=1:operation.nVarOut
    self.operationVariableNameMap(operation.outputVariables(iVar).name) = [];
end
self.operationNameMap{operation.name} = [];
self.clearVariableCache();
end