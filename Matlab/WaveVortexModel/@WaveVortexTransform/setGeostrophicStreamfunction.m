function setGeostrophicStreamfunction(self,psi)
% set a geostrophic streamfunction
%
% Clears A0 by setting a geostrophic streamfunction
% - Topic: Initial conditions — Geostrophic Motions
% - Declaration: setGeostrophicStreamfunction(psi)
% - Parameter psi: function handle that takes three arguments, psi(X,Y,Z)
[X,Y,Z] = ndgrid(self.x,self.y,self.z);
self.A0 = self.transformFromSpatialDomainWithF( (self.f/self.g)*psi(X,Y,Z) );
end