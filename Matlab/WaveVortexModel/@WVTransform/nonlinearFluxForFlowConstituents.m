function [Fp,Fm,F0] = nonlinearFluxForFlowConstituents(self,Uconstituent,gradUconstituent)
% Apply operator T_\omega---defined in (C2) in the manuscript
phase = exp(self.iOmega*(self.t-self.t0));
Apt = self.Ap .* phase;
Amt = self.Am .* conj(phase);
A0t = self.A0;

[ApmUMask,A0UMask] = self.masksForFlowConstituents(Uconstituent);

% Apply operator S---defined in (C4) in the manuscript
Ubar = ApmUMask.*(self.UAp.*Apt + self.UAm.*Amt) + A0UMask.*self.UA0.*A0t;
Vbar = ApmUMask.*(self.VAp.*Apt + self.VAm.*Amt) + A0UMask.*self.VA0.*A0t;
Wbar = ApmUMask.*(self.WAp.*Apt + self.WAm.*Amt);
Nbar = ApmUMask.*(self.NAp.*Apt + self.NAm.*Amt) + A0UMask.*self.NA0.*A0t;

% Finishing applying S, but also compute derivatives at the
% same time
U = self.transformToSpatialDomainWithF(Ubar);
V = self.transformToSpatialDomainWithF(Vbar);
W = self.transformToSpatialDomainWithG(Wbar);
ETA = self.transformToSpatialDomainWithG(Nbar);

% Finishing applying S, but also compute derivatives at the
% same time
[ApmUxMask,A0UxMask] = self.masksForFlowConstituents(gradUconstituent);
Uxbar = ApmUxMask.*(self.UAp.*Apt + self.UAm.*Amt) + A0UxMask.*self.UA0.*A0t;
Vxbar = ApmUxMask.*(self.VAp.*Apt + self.VAm.*Amt) + A0UxMask.*self.VA0.*A0t;
Nxbar = ApmUxMask.*(self.NAp.*Apt + self.NAm.*Amt) + A0UxMask.*self.NA0.*A0t;
[~,Ux,Uy,Uz] = self.transformToSpatialDomainWithFAllDerivatives(Uxbar);
[~,Vx,Vy,Vz] = self.transformToSpatialDomainWithFAllDerivatives(Vxbar);
[~,ETAx,ETAy,ETAz] = self.transformToSpatialDomainWithGAllDerivatives(Nxbar);

% Compute the nonlinear terms in the spatial domain
% (pseudospectral!)
uNL = -U.*Ux - V.*Uy - W.*Uz;
vNL = -U.*Vx - V.*Vy - W.*Vz;
nNL = -U.*ETAx - V.*ETAy - W.*ETAz;
% nNL = -U.*ETAx - V.*ETAy - W.*(ETAz + ETA.*shiftdim(self.dLnN2,-2));


% Now apply the operator S^{-1} and then T_\omega^{-1}
uNLbar = self.transformFromSpatialDomainWithF(uNL);
vNLbar = self.transformFromSpatialDomainWithF(vNL);
nNLbar = self.transformFromSpatialDomainWithG(nNL);

Fp = (self.ApU.*uNLbar + self.ApV.*vNLbar + self.ApN.*nNLbar) .* conj(phase);
Fm = (self.AmU.*uNLbar + self.AmV.*vNLbar + self.AmN.*nNLbar) .* phase;
F0 = self.A0U.*uNLbar + self.A0V.*vNLbar + self.A0N.*nNLbar;
end