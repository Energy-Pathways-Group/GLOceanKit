function varargout = performOperationWithName(self,opName,varargin)
% computes (runs) the operation
%
% - Topic: Internal
modelOp = self.operationNameMap(opName);
varargout = cell(1,modelOp.nVarOut);
[varargout{:}] = self.performOperation(modelOp,varargin{:});
end