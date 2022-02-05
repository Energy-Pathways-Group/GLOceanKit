function [Fp,Fm,F0] = NonlinearFluxForFlowConstituentsAtTime(self,t,Ap,Am,A0,Uconstituent,gradUconstituent)
% Apply operator T_\omega---defined in (C2) in the manuscript
phase = exp(self.iOmega*(t-self.t0));
Ap = Ap .* phase;
Am = Am .* conj(phase);

[ApmUMask,A0UMask] = self.MasksForFlowContinuents(Uconstituent);

% Apply operator S---defined in (C4) in the manuscript
Ubar = ApmUMask.*(self.UAp.*Ap + self.UAm.*Am) + A0UMask.*self.UA0.*A0;
Vbar = ApmUMask.*(self.VAp.*Ap + self.VAm.*Am) + A0UMask.*self.VA0.*A0;
Wbar = ApmUMask.*(self.WAp.*Ap + self.WAm.*Am);
Nbar = ApmUMask.*(self.NAp.*Ap + self.NAm.*Am) + A0UMask.*self.NA0.*A0;

% Finishing applying S, but also compute derivatives at the
% same time
U = self.TransformToSpatialDomainWithF(Ubar);
V = self.TransformToSpatialDomainWithF(Vbar);
W = self.TransformToSpatialDomainWithG(Wbar);
ETA = self.TransformToSpatialDomainWithG(Nbar);

% Finishing applying S, but also compute derivatives at the
% same time
[ApmUxMask,A0UxMask] = self.MasksForFlowContinuents(gradUconstituent);
Uxbar = ApmUxMask.*(self.UAp.*Ap + self.UAm.*Am) + A0UxMask.*self.UA0.*A0;
Vxbar = ApmUxMask.*(self.VAp.*Ap + self.VAm.*Am) + A0UxMask.*self.VA0.*A0;
Nxbar = ApmUxMask.*(self.NAp.*Ap + self.NAm.*Am) + A0UxMask.*self.NA0.*A0;
[~,Ux,Uy,Uz] = self.TransformToSpatialDomainWithFAllDerivatives(Uxbar);
[~,Vx,Vy,Vz] = self.TransformToSpatialDomainWithFAllDerivatives(Vxbar);
[~,ETAx,ETAy,ETAz] = self.TransformToSpatialDomainWithGAllDerivatives(Nxbar);

% Compute the nonlinear terms in the spatial domain
% (pseudospectral!)
uNL = -U.*Ux - V.*Uy - W.*Uz;
vNL = -U.*Vx - V.*Vy - W.*Vz;
nNL = -U.*ETAx - V.*ETAy - W.*(ETAz + ETA.*shiftdim(self.dLnN2,-2));

% Now apply the operator S^{-1} and then T_\omega^{-1}
uNLbar = self.TransformFromSpatialDomainWithF(uNL);
vNLbar = self.TransformFromSpatialDomainWithF(vNL);
nNLbar = self.TransformFromSpatialDomainWithG(nNL);

Fp = (self.ApU.*uNLbar + self.ApV.*vNLbar + self.ApN.*nNLbar) .* conj(phase);
Fm = (self.AmU.*uNLbar + self.AmV.*vNLbar + self.AmN.*nNLbar) .* phase;
F0 = self.A0U.*uNLbar + self.A0V.*vNLbar + self.A0N.*nNLbar;
end