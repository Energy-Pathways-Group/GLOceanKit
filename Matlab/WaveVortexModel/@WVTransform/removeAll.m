function removeAll(self)
% removes all energy from the model
%
% Simply sets Ap, Am, and A0 to zero.
% - Topic: Initial conditions — Waves
self.A0 = 0*self.A0;
self.Ap = 0*self.Ap;
self.Am = 0*self.Am;
end