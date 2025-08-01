function addUVEta(self,U,V,N)
% add $$(u,v,\eta)$$ to the existing values
%
% Add $$(u,v,\eta)$$ to the existing Ap,Am,A0
% - Topic: Initial conditions
% - Declaration: addUVEta(U,V,N)
% - Parameter u: x-component of the fluid velocity
% - Parameter v: y-component of the fluid velocity
% - Parameter n: scaled density anomaly
arguments
    self WVTransform {mustBeNonempty}
    U (:,:,:) double {mustBeNonempty,mustBeReal}
    V (:,:,:) double {mustBeNonempty,mustBeReal}
    N (:,:,:) double {mustBeNonempty,mustBeReal}
end

[Ap,Am,A0] = self.transformUVEtaToWaveVortex(U,V,N);

self.Ap = self.Ap + Ap;
self.Am = self.Am + Am;
self.A0 = self.A0  + A0;

end