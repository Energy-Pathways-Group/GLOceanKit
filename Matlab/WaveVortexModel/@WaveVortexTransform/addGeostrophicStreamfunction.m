function addGeostrophicStreamfunction(self,psi)
[X,Y,Z] = ndgrid(self.x,self.y,self.z);
self.A0 = self.A0 + self.transformFromSpatialDomainWithF( (self.f0/self.g)*psi(X,Y,Z) );
end