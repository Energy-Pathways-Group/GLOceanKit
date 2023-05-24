function setGeostrophicStreamfunction(self,psi)
% set a geostrophic streamfunction
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
% Clears A0 by setting a geostrophic streamfunction
% - Topic: Initial conditions — Geostrophic Motions
% - Declaration: setGeostrophicStreamfunction(psi)
% - Parameter psi: function handle that takes three arguments, psi(X,Y,Z)
[X,Y,Z] = ndgrid(self.x,self.y,self.z);
self.A0 = self.transformFromSpatialDomainWithF( (self.f/self.g)*psi(X,Y,Z) );
end