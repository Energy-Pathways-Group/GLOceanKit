function [Fp,Fm,F0,u_f,v_f,w_f,u_d,v_d] = NonlinearFluxWithFloatsAndDriftersAtTime(self,t,Ap,Am,A0,x_f,y_f,z_f,x_d,y_d,z_d)
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
[U,Ux,Uy,Uz] = self.transformToSpatialDomainWithFAllDerivatives(Ubar);
[V,Vx,Vy,Vz] = self.transformToSpatialDomainWithFAllDerivatives(Vbar);
W = self.transformToSpatialDomainWithG(Wbar);
[~,ETAx,ETAy,ETAz] = self.transformToSpatialDomainWithGAllDerivatives(Nbar);

% Compute the nonlinear terms in the spatial domain
% (pseudospectral!)
uNL = -U.*Ux - V.*Uy - W.*Uz;
vNL = -U.*Vx - V.*Vy - W.*Vz;
nNL = -U.*ETAx - V.*ETAy - W.*ETAz;

% Now apply the operator S^{-1} and then T_\omega^{-1}
uNLbar = self.transformFromSpatialDomainWithF(uNL);
vNLbar = self.transformFromSpatialDomainWithF(vNL);
nNLbar = self.transformFromSpatialDomainWithG(nNL);

Fp = (self.ApU.*uNLbar + self.ApV.*vNLbar + self.ApN.*nNLbar) .* conj(phase);
Fm = (self.AmU.*uNLbar + self.AmV.*vNLbar + self.AmN.*nNLbar) .* phase;
F0 = self.A0U.*uNLbar + self.A0V.*vNLbar + self.A0N.*nNLbar;

[u_f,v_f,w_f] = self.InterpolatedFieldAtPosition(x_f,y_f,z_f,'spline',U,V,W);
[u_d,v_d] = self.InterpolatedFieldAtPosition(x_d,y_d,z_d,'spline',U,V);
end