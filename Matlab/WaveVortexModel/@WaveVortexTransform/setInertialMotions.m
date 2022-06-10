function setInertialMotions(self,u,v)
[~,~,Z] = ndgrid(self.x,self.y,self.z);

Ubar = self.transformFromSpatialDomainWithF( u(Z) );
Vbar = self.transformFromSpatialDomainWithF( v(Z) );
Ap_ = self.ApU.*Ubar + self.ApV.*Vbar;
Am_ = self.AmU.*Ubar + self.AmV.*Vbar;
self.Ap(1,1,:) = Ap_(1,1,:);
self.Am(1,1,:) = Am_(1,1,:);
end