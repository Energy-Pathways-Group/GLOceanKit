function setInertialMotions(self,u,v)
% set inertial motions
%
% Overwrites existing inertial motions with the new values
% - Topic: Initial conditions — Inertial Oscillations
% - Declaration: setInertialMotions(psi)
% - Parameter u: function handle that takes a single argument, u(Z)
% - Parameter v: function handle that takes a single argument, v(Z)
[~,~,Z] = ndgrid(self.x,self.y,self.z);

Ubar = self.transformFromSpatialDomainWithF( u(Z) );
Vbar = self.transformFromSpatialDomainWithF( v(Z) );
Ap_ = self.ApU.*Ubar + self.ApV.*Vbar;
Am_ = self.AmU.*Ubar + self.AmV.*Vbar;
self.Ap(1,1,:) = Ap_(1,1,:);
self.Am(1,1,:) = Am_(1,1,:);
end