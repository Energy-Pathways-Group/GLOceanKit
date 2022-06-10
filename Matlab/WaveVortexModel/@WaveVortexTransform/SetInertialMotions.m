function SetInertialMotions(self,u,v)
[~,~,Z] = ndgrid(self.x,self.y,self.z);

Ubar = self.transformFromSpatialDomainWithF( u(Z) );
Vbar = self.transformFromSpatialDomainWithF( v(Z) );
self.Ap = self.ApU.*Ubar + self.ApV.*Vbar;
self.Am = self.AmU.*Ubar + self.AmV.*Vbar;
warning('SetInertialMotions kills the internal waves')
end