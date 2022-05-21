function [u,v,w] = VelocityFieldAtTime(self, t)
% Return the velocity field, which is the sum of the gridded
% and external/free waves at time t. Note that if you do not
% need w, don't request it and it won't be computed.
if nargout == 3
    [u,v,w] = self.VariableFieldsAtTime(t,'u','v','w');
else
    [u,v] = self.VariableFieldsAtTime(t,'u','v');
end
end