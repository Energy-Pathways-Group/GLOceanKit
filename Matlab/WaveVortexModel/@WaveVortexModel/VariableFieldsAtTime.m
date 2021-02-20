function [varargout] = VariableFieldsAtTime(self, t, varargin)
% Primary method for accessing the dynamical variables on the
% internal grid.
%
% Valid variable options are 'u', 'v', 'w', 'rho_prime', and
% 'zeta'.
varargout = cell(size(varargin));
[varargout{:}] = self.InternalVariableFieldsAtTime(t, varargin{:});
end