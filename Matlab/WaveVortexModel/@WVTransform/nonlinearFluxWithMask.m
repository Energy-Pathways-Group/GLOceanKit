function [Fp,Fm,F0] = nonlinearFluxWithMask(self,mask)
% returns the flux of each coefficient as determined by the nonlinear flux
% operation and the given mask.
%
% The mask is applied to the coefficients Ap,Am,A0 before computing the
% nonlinear flux. This is useful for zeroing wavenumbers at given total
% wavenumber or frequency, for example.
%
% The nonlinear flux used is the unforced, invicid equations.
%
% - Topic: Nonlinear flux and energy transfers
% - Declaration: [Fp,Fm,F0] = nonlinearFluxWithMask(mask)
% - Parameter mask: mask applied to all constituents
% - Returns Fp: flux into the Ap coefficients
% - Returns Fm: flux into the Am coefficients
% - Returns F0: flux into the A0 coefficients
arguments
    self WVTransform {mustBeNonempty}
    mask (:,:) double {mustBeNonempty,mustBeReal}
end

% Apply operator T_\omega---defined in (C2) in the manuscript
phase = exp(self.iOmega*(self.t-self.t0));
Apt = mask .* self.Ap .* phase;
Amt = mask .* self.Am .* conj(phase);
A0t = mask .* self.A0;

% Apply operator S---defined in (C4) in the manuscript
Ubar = self.UAp.*Apt + self.UAm.*Amt + self.UA0.*A0t;
Vbar = self.VAp.*Apt + self.VAm.*Amt + self.VA0.*A0t;
Wbar = self.WAp.*Apt + self.WAm.*Amt;
Nbar = self.NAp.*Apt + self.NAm.*Amt + self.NA0.*A0t;

% Finishing applying S, but also compute derivatives at the same time
[U,Ux,Uy,Uz] = self.transformToSpatialDomainWithFAllDerivatives(Ubar);
[V,Vx,Vy,Vz] = self.transformToSpatialDomainWithFAllDerivatives(Vbar);
W = self.transformToSpatialDomainWithG(Wbar);
[ETA,ETAx,ETAy,ETAz] = self.transformToSpatialDomainWithGAllDerivatives(Nbar);

% Compute the nonlinear terms in the spatial domain (pseudospectral!)
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

[ApNL,AmNL,A0NL] = self.transformUVEtaToWaveVortex(uNL,vNL,nNL);

Fp = ApNL .* conj(phase);
Fm = AmNL .* phase;
F0 = A0NL;
end