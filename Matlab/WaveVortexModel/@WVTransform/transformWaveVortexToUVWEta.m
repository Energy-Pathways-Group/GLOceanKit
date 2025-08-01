function [U,V,W,N] = transformWaveVortexToUVWEta(self,Ap,Am,A0,t)
% transform wave-vortex coefficients $$(A_+,A_-,A_0)$$ to fluid variables $$(u,v,\eta)$$.
%
% This function is the inverse WVTransform. It is a
% [linear
% transformation](/mathematical-introduction/transformations.html)
% denoted $$\mathcal{L}$$.
%
% This function is not intended to be used directly (although
% you can), and is kept here to demonstrate a simple
% implementation of the transformation. Instead, you should
% initialize the WVTransform using one of the
% initialization functions.
%
% - Topic: Operations — Transformations
% - Declaration: [u,v,w,n] = transformWaveVortexToUVWEta(self,Ap,Am,A0,t)
% - Parameter Ap: positive wave coefficients at reference time t0
% - Parameter Am: negative wave coefficients at reference time t0
% - Parameter A0: geostrophic coefficients at reference time t0
% - Parameter t: (optional) time of observations
% - Returns u: x-component of the fluid velocity
% - Returns v: y-component of the fluid velocity
% - Returns w: z-component of the fluid velocity
% - Returns n: scaled density anomaly

if nargin == 5
    phase = exp(self.iOmega*(t-self.t0));
    Ap = Ap .* phase;
    Am = Am .* conj(phase);
end

% This is the 'S' operator (C4) in the manuscript
U = self.transformToSpatialDomainWithF(Apm=self.UAp.*Ap + self.UAm.*Am, A0=self.UA0.*A0);
V = self.transformToSpatialDomainWithF(Apm=self.VAp.*Ap + self.VAm.*Am, A0=self.VA0.*A0);
W = self.transformToSpatialDomainWithG(Apm=self.WAp.*Ap + self.WAm.*Am);
N = self.transformToSpatialDomainWithG(Apm=self.NAp.*Ap + self.NAm.*Am, A0=self.NA0.*A0);
end