function varargout = performOperationWithName(self,opName)
% computes (runs) the operation
%
% - Topic: Internal
arguments (Input)
    self WVTransform
    opName char
end
arguments (Output,Repeating)
    varargout
end
modelOp = self.operationNameMap{opName};
varargout = cell(1,modelOp.nVarOut);
[varargout{:}] = self.performOperation(modelOp);
end