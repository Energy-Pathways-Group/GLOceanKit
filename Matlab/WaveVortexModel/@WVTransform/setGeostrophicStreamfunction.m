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
% N^2 \eta = \frac{g}{\rho_0} \rho = - f \frac{\partial \psi}{\partial z}
% $$
%
% Clears A0 by setting a geostrophic streamfunction
% - Topic: Initial conditions — Geostrophic Motions
% - Declaration: setGeostrophicStreamfunction(psi)
% - Parameter psi: function handle that takes three arguments, psi(X,Y,Z)
[X,Y,Z] = self.xyzGrid;
psi_bar = self.transformFromSpatialDomainWithFourier((self.f/self.g)*psi(X,Y,Z) );
% psi_barz = self.diffZF(psi_bar);
% psi_bar(1,1,:) = 0;
self.A0 = self.transformFromSpatialDomainWithFg( psi_bar);

% a = -psi_barz(1,1,:)./shiftdim(self.N2,-2);
% psi_bar0z = self.transformFromSpatialDomainWithGmda(a);
% self.A0(1,1,:) = psi_bar0z;

end