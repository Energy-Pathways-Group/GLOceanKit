function [Ap,Am,A0] = transformUVEtaToWaveVortex(self,U,V,N)
% transform fluid variables $$(u,v,\eta)$$ to wave-vortex coefficients $$(A_+,A_-,A_0)$$.
%
% This function **is** the WVTransform. It is a [linear
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
% - Declaration: [Ap,Am,A0] = transformUVEtaToWaveVortex(U,V,N)
% - Parameter u: x-component of the fluid velocity
% - Parameter v: y-component of the fluid velocity
% - Parameter n: scaled density anomaly
% - Returns Ap: positive wave coefficients at reference time t0
% - Returns Am: negative wave coefficients at reference time t0
% - Returns A0: geostrophic coefficients at reference time t0
u_hat = self.transformFromSpatialDomainWithFourier(U);
v_hat = self.transformFromSpatialDomainWithFourier(V);
n_hat = self.transformFromSpatialDomainWithFourier(N);

iK = sqrt(-1)*repmat(shiftdim(self.k,-1),self.Nz,1);
iL = sqrt(-1)*repmat(shiftdim(self.l,-1),self.Nz,1);

n_bar = self.transformFromSpatialDomainWithGg(n_hat);
zeta_bar = self.transformFromSpatialDomainWithFg(iK .* v_hat - iL .* u_hat);
A0 = self.A0Z.*zeta_bar + self.A0N.*n_bar;

delta_bar = self.transformWithG_wg(self.h_0.*self.transformFromSpatialDomainWithFg(iK .* u_hat + iL .* v_hat));
nw_bar = self.transformWithG_wg(n_bar - self.NA0.*A0);
Ap = self.ApmD .* delta_bar + self.ApmN .* nw_bar;
Am = self.ApmD .* delta_bar - self.ApmN .* nw_bar;

Ap(:,1) = self.transformFromSpatialDomainWithFio(u_hat(:,1) - sqrt(-1)*v_hat(:,1))/2;
Am(:,1) = conj(Ap(:,1));

if nargin == 5
    phase = exp(-self.iOmega*(t-self.t0));
    Ap = Ap .* phase;
    Am = Am .* conj(phase);
end
end