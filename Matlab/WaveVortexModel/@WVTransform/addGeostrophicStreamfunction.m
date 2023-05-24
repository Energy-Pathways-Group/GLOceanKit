function addGeostrophicStreamfunction(self,psi)
% add a geostrophic streamfunction to existing geostrophic motions
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
% The geostrophic streamfunction is added to the existing values in `A0`
% - Topic: Initial conditions — Geostrophic Motions
% - Declaration: addGeostrophicStreamfunction(psi)
% - Parameter psi: function handle that takes three arguments, psi(X,Y,Z)
[X,Y,Z] = ndgrid(self.x,self.y,self.z);
self.A0 = self.A0 + self.transformFromSpatialDomainWithF( (self.f/self.g)*psi(X,Y,Z) );
end