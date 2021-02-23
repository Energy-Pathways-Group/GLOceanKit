function [Fp,Fm,F0,u,v,w] = NonlinearFluxWithParticlesAtTime(self,t,Ap,Am,A0,x,y,z)
% Apply operator T_\omega---defined in (C2) in the manuscript
phase = exp(self.iOmega*(t-self.t0));
Ap = Ap .* phase;
Am = Am .* conj(phase);

% Apply operator S---defined in (C4) in the manuscript
Ubar = self.UAp.*Ap + self.UAm.*Am + self.UA0.*A0;
Vbar = self.VAp.*Ap + self.VAm.*Am + self.VA0.*A0;
Wbar = self.WAp.*Ap + self.WAm.*Am;
Nbar = self.NAp.*Ap + self.NAm.*Am + self.NA0.*A0;

% Finishing applying S, but also compute derivatives at the
% same time
[U,Ux,Uy,Uz] = self.TransformToSpatialDomainWithFAllDerivatives(Ubar);
[V,Vx,Vy,Vz] = self.TransformToSpatialDomainWithFAllDerivatives(Vbar);
W = self.TransformToSpatialDomainWithG(Wbar);
[~,ETAx,ETAy,ETAz] = self.TransformToSpatialDomainWithGAllDerivatives(Nbar);

% Compute the nonlinear terms in the spatial domain
% (pseudospectral!)
uNL = -U.*Ux - V.*Uy - W.*Uz;
vNL = -U.*Vx - V.*Vy - W.*Vz;
nNL = -U.*ETAx - V.*ETAy - W.*ETAz;

% Now apply the operator S^{-1} and then T_\omega^{-1}
uNLbar = self.TransformFromSpatialDomainWithF(uNL);
vNLbar = self.TransformFromSpatialDomainWithF(vNL);
nNLbar = self.TransformFromSpatialDomainWithG(nNL);

Fp = (self.ApU.*uNLbar + self.ApV.*vNLbar + self.ApN.*nNLbar) .* conj(phase);
Fm = (self.AmU.*uNLbar + self.AmV.*vNLbar + self.AmN.*nNLbar) .* phase;
F0 = self.A0U.*uNLbar + self.A0V.*vNLbar + self.A0N.*nNLbar;

[u,v,w] = self.InterpolatedFieldAtPosition(x,y,z,'spline',U,V,W);
end