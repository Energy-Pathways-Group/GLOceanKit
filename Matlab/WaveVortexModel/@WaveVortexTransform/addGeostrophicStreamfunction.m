function addGeostrophicStreamfunction(self,psi)
[X,Y,Z] = ndgrid(self.x,self.y,self.z);
A0_ = (1./self.PP) .* self.transformFromSpatialDomainWithF( (self.f0/self.g)*psi(X,Y,Z) );
self.A0 = self.A0 + A0_;
end