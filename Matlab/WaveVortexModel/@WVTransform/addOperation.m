function addOperation(self,operation,options)
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
    operation (1,:) WVOperation {mustBeNonempty}
    options.shouldOverwriteExisting logical = false
    options.shouldSuppressWarning logical = false
end
for iOp=1:length(operation)
    isExisting = 0;
    for iVar=1:length(operation(iOp).outputVariables)
        isKnownDimension = ismember(operation(iOp).outputVariables(iVar).dimensions,self.annotatedDimensionNames);
        if any(~isKnownDimension)
            error("Unable to find dimension (%s) for the numeric property %s",strjoin(string(operation(iOp).outputVariables(iVar).dimensions(~isKnownDimension)),","),operation(iOp).outputVariables(iVar).name);
        end

        isKnownProperty = ismember(operation(iOp).outputVariables(iVar).name,self.annotatedPropertyNames);
        if any(isKnownProperty)
            existingVar = self.propertyAnnotationWithName(operation(iOp).outputVariables(iVar).name);
            if isa(existingVar,'WVVariableAnnotation')
                isExisting = 1;
            end
        end
    end

    if isExisting == 1
        message1 = strcat(existingVar.modelOp.name,' with variables {');
        for jVar=1:existingVar.modelOp.nVarOut
            message1 = strcat(message1,existingVar.modelOp.outputVariables(jVar).name,',');
        end
        message1 = strcat(message1,'}');

        message2 = strcat(operation(iOp).name,' with variables {');
        for jVar=1:operation(iOp).nVarOut
            message2 = strcat(message2,operation(iOp).outputVariables(jVar).name,',');
        end
        message2 = strcat(message2,'}');
        if options.shouldOverwriteExisting == 0
            error('A variable with the same name already exists! You attempted to replace the operation %s with the operation %s. If you are sure you want to do this, call wvt.addOperation(newOp,overwriteExisting=1).', message1,message2);
        else
            self.removeOperation(existingVar.modelOp);
            if options.shouldSuppressWarning == false
                fprintf('The operation %s has been removed and the operation %s has been added.\n',message1,message2);
            end
        end
    end

    % Now go ahead and actually add the operation and its variables
    self.addPropertyAnnotation(operation(iOp).outputVariables);
    for iVar=1:length(operation(iOp).outputVariables)
        self.operationVariableNameMap(operation(iOp).outputVariables(iVar).name) = operation(iOp).outputVariables(iVar);
    end
    self.operationNameMap{operation(iOp).name} = operation(iOp);

end
end