function varargout = performOperationWithName(self,opName,varargin)
% computes (runs) the operation
%
% - Topic: Internal
modelOp = self.operationNameMap(opName);
varargout = cell(1,modelOp{1}.nVarOut);
[varargout{:}] = self.performOperation(modelOp{1},varargin{:});
end