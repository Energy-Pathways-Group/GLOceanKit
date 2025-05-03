function initWithUVRho(self,U,V,RHO)
% initialize with fluid variables $$(u,v,\rho)$$
%
% Replaces the variables Ap,Am,A0 with those computed from $$(u,v,\rho_e)$$.
% - Topic: Initial conditions
% - Declaration: initWithUVRho(U,V,RHO)
% - Parameter u: x-component of the fluid velocity
% - Parameter v: y-component of the fluid velocity
% - Parameter rho: density anomaly

arguments
    self WVTransform {mustBeNonempty}
    U (:,:,:) double {mustBeNonempty,mustBeReal}
    V (:,:,:) double {mustBeNonempty,mustBeReal}
    RHO (:,:,:) double {mustBeNonempty,mustBeReal}
end

[self.Ap,self.Am,self.A0] = self.transformUVEtaToWaveVortex(U,V,(self.g/self.rho0)*RHO./shiftdim(self.N2,-2));

end