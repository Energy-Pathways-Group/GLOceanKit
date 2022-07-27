function initWithInertialMotions(self,u,v)
% initialize with inertial motions
%
% Clears variables Ap,Am,A0 and then sets inertial motions
% - Topic: Initial conditions — Inertial Oscillations
% - Declaration: initWithInertialMotions(u,v)
% - Parameter u: function handle that takes a single argument, u(Z)
% - Parameter v: function handle that takes a single argument, v(Z)
self.Ap = zeros(size(self.Ap));
self.Am = zeros(size(self.Am));
self.A0 = zeros(size(self.A0));
self.setInertialMotions(u,v);
end