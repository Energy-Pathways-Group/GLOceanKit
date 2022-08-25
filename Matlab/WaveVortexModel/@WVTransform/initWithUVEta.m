function initWithUVEta(self,U,V,N,options)
% initialize with fluid variables $$(u,v,\eta)$$
%
% Replaces the variables Ap,Am,A0 with those computed from $$(u,v,\eta)$$.
% If a time t is specified, the wvt is set to that time.
% - Topic: Initial conditions
% - Declaration: initWithUVEta(U,V,N,t)
% - Parameter u: x-component of the fluid velocity
% - Parameter v: y-component of the fluid velocity
% - Parameter n: scaled density anomaly
% - Parameter t: (optional) time of observations

arguments
    self WVTransform {mustBeNonempty}
    U (:,:,:) double {mustBeNonempty,mustBeReal}
    V (:,:,:) double {mustBeNonempty,mustBeReal}
    N (:,:,:) double {mustBeNonempty,mustBeReal}
    options.t (1,1) double = 0
end

if isfield(options,'t')
    self.t = options.t;
    [Ap,Am,A0] = self.transformUVEtaToWaveVortex(U,V,N,options.t);
else
    [Ap,Am,A0] = self.transformUVEtaToWaveVortex(U,V,N,self.t);
end

self.Ap = Ap;
self.Am = Am;
self.A0 = A0;

end