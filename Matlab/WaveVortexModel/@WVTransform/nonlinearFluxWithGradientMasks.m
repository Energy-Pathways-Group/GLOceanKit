function [Fp,Fm,F0] = nonlinearFluxWithGradientMasks(self,ApUMask,AmUMask,A0UMask,ApUxMask,AmUxMask,A0UxMask)
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
    ApUMask (:,:) double {mustBeNonempty,mustBeReal}
    AmUMask (:,:) double {mustBeNonempty,mustBeReal}
    A0UMask (:,:) double {mustBeNonempty,mustBeReal}
    ApUxMask (:,:) double {mustBeNonempty,mustBeReal}
    AmUxMask (:,:) double {mustBeNonempty,mustBeReal}
    A0UxMask (:,:) double {mustBeNonempty,mustBeReal}
end

% Apply operator T_\omega---defined in (C2) in the manuscript
Apt = self.Apt;
Amt = self.Amt;
A0t = self.A0t;

% Apply operator S---defined in (C4) in the manuscript
% Finishing applying S, but also compute derivatives at the
% same time
U = self.transformToSpatialDomainWithF(Apm=(ApUMask.*self.UAp.*Apt + AmUMask.*self.UAm.*Amt), A0=A0UMask.*self.UA0.*A0t);
V = self.transformToSpatialDomainWithF(Apm=(ApUMask.*self.VAp.*Apt + AmUMask.*self.VAm.*Amt), A0=A0UMask.*self.VA0.*A0t);
W = self.transformToSpatialDomainWithG(Apm=(ApUMask.*self.WAp.*Apt + AmUMask.*self.WAm.*Amt));

% Finishing applying S, but also compute derivatives at the
% same time
[~,Ux,Uy,Uz] = self.transformToSpatialDomainWithFAllDerivatives(        Apm=(ApUxMask.*self.UAp.*Apt + AmUxMask.*self.UAm.*Amt),A0=A0UxMask.*self.UA0.*A0t);
[~,Vx,Vy,Vz] = self.transformToSpatialDomainWithFAllDerivatives(        Apm=(ApUxMask.*self.VAp.*Apt + AmUxMask.*self.VAm.*Amt),A0=A0UxMask.*self.VA0.*A0t);
[ETA,ETAx,ETAy,ETAz] = self.transformToSpatialDomainWithGAllDerivatives(Apm=(ApUxMask.*self.NAp.*Apt + AmUxMask.*self.NAm.*Amt),A0=A0UxMask.*self.NA0.*A0t);

% Compute the nonlinear terms in the spatial domain
% (pseudospectral!)
uNL = -U.*Ux - V.*Uy - W.*Uz;
vNL = -U.*Vx - V.*Vy - W.*Vz;
if isa(self,'WVTransformConstantStratification')
    nNL = -U.*ETAx - V.*ETAy - W.*ETAz;
    [Fp,Fm,F0] = self.transformUVEtaToWaveVortex(uNL,vNL,nNL);
elseif isa(self,'WVTransformHydrostatic')
    nNL = -U.*ETAx - V.*ETAy - W.*(ETAz + ETA.*shiftdim(self.dLnN2,-2));
    [Fp,Fm,F0] = self.transformUVEtaToWaveVortex(uNL,vNL,nNL);
elseif isa(self,'WVTransformBoussinesq')
    [~,Wx,Wy,Wz] = self.transformToSpatialDomainWithGAllDerivatives(Apm=(ApUxMask.*self.WAp.*Apt + AmUxMask.*self.WAm.*Amt));
    wNL = -U.*Wx - V.*Wy - W.*Wz;
    nNL = -U.*ETAx - V.*ETAy - W.*(ETAz + ETA.*shiftdim(self.dLnN2,-2));
    [Fp,Fm,F0] = self.transformUVWEtaToWaveVortex(uNL,vNL,wNL,nNL);
else
    error('WVTransform not recognized.')
end

end