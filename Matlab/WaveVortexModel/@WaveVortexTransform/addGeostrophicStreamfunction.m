function addGeostrophicStreamfunction(self,psi)
% add a geostrophic streamfunction to existing geostrophic motions
%
% The geostrophic streamfunction is added to the existing values in `A0`
% - Topic: Initial conditions — Geostrophic Motions
% - Declaration: addGeostrophicStreamfunction(psi)
% - Parameter psi: function handle that takes three arguments, psi(X,Y,Z)
[X,Y,Z] = ndgrid(self.x,self.y,self.z);
self.A0 = self.A0 + self.transformFromSpatialDomainWithF( (self.f/self.g)*psi(X,Y,Z) );
end