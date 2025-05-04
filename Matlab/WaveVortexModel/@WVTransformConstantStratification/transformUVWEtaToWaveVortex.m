function [Ap,Am,A0] = transformUVWEtaToWaveVortex(self,U,V,W,N)
% transform momentum variables $$(u,v,w,\eta)$$ to wave-vortex coefficients $$(A_+,A_-,A_0)$$.
%
% This function tuned for constant stratification.
%
% - Topic: Operations — Transformations
% - Declaration: [Ap,Am,A0] = transformUVWEtaToWaveVortex(U,V,N)
% - Parameter u: x-component of the momentum
% - Parameter v: y-component of the momentum
% - Parameter w: y-component of the momentum
% - Parameter n: scaled density anomaly
% - Parameter t: (optional) time of observations
% - Returns Ap: positive wave coefficients at reference time t0
% - Returns Am: negative wave coefficients at reference time t0
% - Returns A0: geostrophic coefficients at reference time t0
u_hat = self.transformFromSpatialDomainWithFourier(U);
v_hat = self.transformFromSpatialDomainWithFourier(V);
w_hat = self.transformFromSpatialDomainWithFourier(W);
n_hat = self.transformFromSpatialDomainWithFourier(N);

iK = sqrt(-1)*repmat(shiftdim(self.k,-1),self.Nz,1);
iL = sqrt(-1)*repmat(shiftdim(self.l,-1),self.Nz,1);

n_bar = self.transformFromSpatialDomainWithGg(n_hat);
zeta_bar = self.transformFromSpatialDomainWithFg(iK .* v_hat - iL .* u_hat);
A0 = self.A0Z.*zeta_bar + self.A0N.*n_bar;
nw_bar = self.transformWithG_wg(n_bar - self.NA0.*A0);

delta_bar = self.ApmD_scaled .* (self.DCT * (self.cos_alpha .* u_hat + self.sin_alpha .* v_hat));
w_bar = self.ApmW_scaled .* (self.DST * w_hat);
Ap = delta_bar + w_bar + self.ApmN .* nw_bar;
Am = delta_bar + w_bar - self.ApmN .* nw_bar;

Ap(:,1) = self.transformFromSpatialDomainWithFio(u_hat(:,1) - sqrt(-1)*v_hat(:,1))/2;
Am(:,1) = conj(Ap(:,1));

Ap = Ap .* self.conjPhase;
Am = Am .* self.phase;
end