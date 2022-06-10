function SetGeostrophicStreamfunction(self,psi)
[X,Y,Z] = ndgrid(self.x,self.y,self.z);
self.A0 = (1./self.PP) .* self.transformFromSpatialDomainWithF( (self.f0/self.g)*psi(X,Y,Z) );
warning('SetGeostrophicStreamfunction is only valid for hydrostatic modes! Still need to generalize.');
end