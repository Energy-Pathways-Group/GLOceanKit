function initWithGeostrophicStreamfunction(self,psi)
% initialize with a geostrophic streamfunction
%
% The geostrophic streamfunction, $$\psi$$, is defined such that
%
% $$
% u= - \frac{\partial \psi}{\partial y}
% $$
% 
% $$
% v=\frac{\partial \psi}{\partial x}
% $$
% 
% $$
% N^2 \eta = - f \frac{\partial \psi}{\partial z}
% $$
%
% Clears variables Ap,Am,A0 and then sets the geostrophic streamfunction
% - Topic: Initial conditions — Geostrophic Motions
% - Declaration: initWithGeostrophicStreamfunction(psi)
% - Parameter psi: function handle that takes three arguments, psi(X,Y,Z)
self.Ap = zeros(size(self.Ap));
self.Am = zeros(size(self.Am));
self.A0 = zeros(size(self.A0));
self.setGeostrophicStreamfunction(psi);
end