function [Fp,Fm,F0] = nonlinearFluxWithGradientMasks(self,ApmUMask,A0UMask,ApmUxMask,A0UxMask)
% returns the flux of each coefficient as determined by the nonlinear flux operation
%
% The masks are applied to the coefficients Ap,Am,A0 before computing the
% nonlinear flux, $$\vec{u} \cdot \nabla \vec{u}$$. This function offers
% more fine-grained control than -nonlinearFluxWithMask.
%
% The nonlinear flux used is the unforced, invicid equations.
%
% - Topic: Nonlinear flux and energy transfers
% - Declaration: [Fp,Fm,F0] = nonlinearFluxWithGradientMasks(ApmUMask,A0UMask,ApmUxMask,A0UxMask)
% - Parameter ApmUMask: mask for the wave portion of $$\vec{u}$$ in $$\vec{u} \cdot \nabla \vec{u}$$
% - Parameter A0UMask: mask for the geostrophic portion of $$\vec{u}$$ in $$\vec{u} \cdot \nabla \vec{u}$$
% - Parameter ApmUxMask: mask for the wave portion of $$\nabla \vec{u}$$ in $$\vec{u} \cdot \nabla \vec{u}$$
% - Parameter A0UxMask: mask for the geostrophic portion of $$\nabla \vec{u}$$ in $$\vec{u} \cdot \nabla \vec{u}$$
% - Returns Fp: flux into the Ap coefficients
% - Returns Fm: flux into the Am coefficients
% - Returns F0: flux into the A0 coefficients
arguments
    self WVTransform {mustBeNonempty}
    ApmUMask (:,:,:) double {mustBeNonempty,mustBeReal}
    A0UMask (:,:,:) double {mustBeNonempty,mustBeReal}
    ApmUxMask (:,:,:) double {mustBeNonempty,mustBeReal}
    A0UxMask (:,:,:) double {mustBeNonempty,mustBeReal}
end

% Apply operator T_\omega---defined in (C2) in the manuscript
phase = exp(self.iOmega*(self.t-self.t0));
Apt = self.Ap .* phase;
Amt = self.Am .* conj(phase);
A0t = self.A0;

% Apply operator S---defined in (C4) in the manuscript
Ubar = ApmUMask.*(self.UAp.*Apt + self.UAm.*Amt) + A0UMask.*self.UA0.*A0t;
Vbar = ApmUMask.*(self.VAp.*Apt + self.VAm.*Amt) + A0UMask.*self.VA0.*A0t;
Wbar = ApmUMask.*(self.WAp.*Apt + self.WAm.*Amt);
% Nbar = ApmUMask.*(self.NAp.*Apt + self.NAm.*Amt) + A0UMask.*self.NA0.*A0t;

% Finishing applying S, but also compute derivatives at the
% same time
U = self.transformToSpatialDomainWithF(Ubar);
V = self.transformToSpatialDomainWithF(Vbar);
W = self.transformToSpatialDomainWithG(Wbar);
% ETA = self.transformToSpatialDomainWithG(Nbar);

% Finishing applying S, but also compute derivatives at the
% same time
Uxbar = ApmUxMask.*(self.UAp.*Apt + self.UAm.*Amt) + A0UxMask.*self.UA0.*A0t;
Vxbar = ApmUxMask.*(self.VAp.*Apt + self.VAm.*Amt) + A0UxMask.*self.VA0.*A0t;
Nxbar = ApmUxMask.*(self.NAp.*Apt + self.NAm.*Amt) + A0UxMask.*self.NA0.*A0t;
[~,Ux,Uy,Uz] = self.transformToSpatialDomainWithFAllDerivatives(Uxbar);
[~,Vx,Vy,Vz] = self.transformToSpatialDomainWithFAllDerivatives(Vxbar);
[ETA,ETAx,ETAy,ETAz] = self.transformToSpatialDomainWithGAllDerivatives(Nxbar);

% Compute the nonlinear terms in the spatial domain
% (pseudospectral!)
uNL = -U.*Ux - V.*Uy - W.*Uz;
vNL = -U.*Vx - V.*Vy - W.*Vz;
if isa(self,'WVTransformConstantStratification')
    nNL = -U.*ETAx - V.*ETAy - W.*ETAz;
elseif isa(self,'WVTransformHydrostatic')
    nNL = -U.*ETAx - V.*ETAy - W.*(ETAz + ETA.*shiftdim(self.dLnN2,-2));
else
    nNL = -U.*ETAx - V.*ETAy - W.*(ETAz + ETA.*shiftdim(self.dLnN2,-2));
    warning('WVTransform not recognized.')
end

% Now apply the operator S^{-1} and then T_\omega^{-1}
uNLbar = self.transformFromSpatialDomainWithF(uNL);
vNLbar = self.transformFromSpatialDomainWithF(vNL);
nNLbar = self.transformFromSpatialDomainWithG(nNL);

Fp = (self.ApU.*uNLbar + self.ApV.*vNLbar + self.ApN.*nNLbar) .* conj(phase);
Fm = (self.AmU.*uNLbar + self.AmV.*vNLbar + self.AmN.*nNLbar) .* phase;
F0 = self.A0U.*uNLbar + self.A0V.*vNLbar + self.A0N.*nNLbar;
end