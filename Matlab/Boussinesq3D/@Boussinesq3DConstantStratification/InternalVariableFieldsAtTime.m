function [varargout] = InternalVariableFieldsAtTime(self, t, varargin)
% This is the primary function for computing the internal
% gridded dynamical variables. It tries to be somewhat
% optimized, by only computing the phase once, and only
% computing the requested variables.
%
% Valid variable options are 'u', 'v', 'w', 'rho_prime', and
% 'eta'.
if length(varargin) < 1
    return;
end

phase = exp(self.iOmega*(t-self.t0));

Ap = self.Ap .* phase;
Am = self.Am .* conj(phase);
A0 = self.A0;

varargout = cell(size(varargin));
for iArg=1:length(varargin)
    if ( strcmp(varargin{iArg}, 'u') )
        varargout{iArg} = self.TransformToSpatialDomainWithF(self.UAp.*Ap + self.UAm.*Am + self.UA0.*A0);
    elseif ( strcmp(varargin{iArg}, 'v') )
        varargout{iArg} = self.TransformToSpatialDomainWithF(self.VAp.*Ap + self.VAm.*Am + self.VA0.*A0);
    elseif ( strcmp(varargin{iArg}, 'w') )
        varargout{iArg} = self.TransformToSpatialDomainWithG(self.WAp.*Ap + self.WAm.*Am);
    elseif ( strcmp(varargin{iArg}, 'p') )
        varargout{iArg} = self.rho0*self.g*self.TransformToSpatialDomainWithF(self.NAp.*Ap + self.NAm.*Am + self.NA0.*A0);
    elseif ( strcmp(varargin{iArg}, 'rho_prime') )
        varargout{iArg} = (self.rho0/9.81)*self.N0*self.N0*self.TransformToSpatialDomainWithG(self.NAp.*Ap + self.NAm.*Am + self.NA0.*A0);
    elseif ( strcmp(varargin{iArg}, 'eta') )
        varargout{iArg} = self.TransformToSpatialDomainWithG(self.NAp.*Ap + self.NAm.*Am + self.NA0.*A0);
    else
        error('Invalid option. You may request u, v, w, rho_prime, or eta.');
    end
end
end