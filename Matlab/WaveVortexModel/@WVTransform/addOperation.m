function addOperation(self,transformOperation,options)
% add a WVOperation
%
% Several things happen when adding an operation.
% 1. We check that dimensions exist for all output variables
% produced by this operation.
% 2. We see if there are any existing output variables with the
% same name.
%   2a. We remove the operation that produced the existing
%   variables, if it exists.
% 3. We map each new variable to this operation variableAnnotationNameMap
% 4. Map each operation name to the operation
%
% In our revision,
%   - The variableAnnotationNameMap will map the name to the
%   variable annotation
%   - The operationVariableNameMap will map the name to the
%   operation
%
% - Topic: Utility function — Metadata
arguments
    self WVTransform {mustBeNonempty}
    transformOperation (1,:) WVOperation {mustBeNonempty}
    options.overwriteExisting double = 0
end
for iOp=1:length(transformOperation)
    isExisting = 0;
    for iVar=1:length(transformOperation(iOp).outputVariables)
        if any(~isKey(self.dimensionAnnotationNameMap,transformOperation(iOp).outputVariables(iVar).dimensions))
            error("Unable to find at least one of the dimensions for variable %s",transformOperation(iOp).outputVariables(iVar).name);
        end

        if isKey(self.variableAnnotationNameMap,transformOperation(iOp).outputVariables(iVar).name)
            isExisting = 1;
            existingVar = self.variableAnnotationWithName(transformOperation(iOp).outputVariables(iVar).name);
        end
    end

    if isExisting == 1
        message1 = strcat(existingVar.modelOp.name,' with variables {');
        for jVar=1:existingVar.modelOp.nVarOut
            message1 = strcat(message1,existingVar.modelOp.outputVariables(jVar).name,',');
        end
        message1 = strcat(message1,'}');

        message2 = strcat(transformOperation(iOp).name,' with variables {');
        for jVar=1:transformOperation(iOp).nVarOut
            message2 = strcat(message2,transformOperation(iOp).outputVariables(jVar).name,',');
        end
        message2 = strcat(message2,'}');
        if options.overwriteExisting == 0
            error('A variable with the same name already exists! You attempted to replace the operation %s with the operation %s. If you are sure you want to do this, call wvt.addOperation(newOp,overwriteExisting=1).', message1,message2);
        else
            self.removeOperation(existingVar.modelOp);
            fprintf('The operation %s has been removed and the operation %s has been added.\n',message1,message2);
        end
    end

    % Now go ahead and actually add the operation and its variables
    self.addVariableAnnotations(transformOperation(iOp).outputVariables);
    for iVar=1:length(transformOperation(iOp).outputVariables)
        self.operationVariableNameMap(transformOperation(iOp).outputVariables(iVar).name) = transformOperation(iOp).outputVariables(iVar);
    end
    self.operationNameMap(transformOperation(iOp).name) = {transformOperation(iOp)};

end
end