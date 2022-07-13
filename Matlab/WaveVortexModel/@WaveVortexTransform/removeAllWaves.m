function removeAllWaves(self)
% removes all wave from the model, including inertial oscillations
%
% Simply sets Ap and Am to zero.
% - Topic: Initial conditions — Waves
self.Ap = 0*self.Ap;
self.Am = 0*self.Am;
end