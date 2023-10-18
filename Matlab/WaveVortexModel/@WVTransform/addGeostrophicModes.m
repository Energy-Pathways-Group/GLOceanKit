function [k,l] = addGeostrophicModes(self, vortexproperties)
% add amplitudes of the given geostrophic modes
%
% Add new amplitudes to any existing amplitudes
% - Topic: Initial conditions — Geostrophic Motions
% - Declaration: [k,l] = addGeostrophicModes(self)
% - Parameter k: integer index, (k0 > -Nx/2 && k0 < Nx/2)
% - Parameter l: integer index, (l0 > -Ny/2 && l0 < Ny/2)
% - Parameter j: integer index, (j0 >= 1 && j0 <= nModes), unless k=l=j=0
% - Parameter phi: phase in radians, (0 <= phi <= 2*pi)
% - Parameter u: fluid velocity u (m/s)
% - Returns k: wavenumber k of the waves (radians/m)
% - Returns l: wavenumber l of the waves (radians/m)
arguments
    self WVTransform {mustBeNonempty}
    vortexproperties.k (:,1) double
    vortexproperties.l (:,1) double
    vortexproperties.j (:,1) double
    vortexproperties.phi (:,1) double
    vortexproperties.u (:,1) double
end
kMode = vortexproperties.k;
lMode = vortexproperties.l;
jMode = vortexproperties.j;
phi = vortexproperties.phi;
u = vortexproperties.u;

[kIndex,lIndex,jIndex,A0Amp] = self.geostrophicCoefficientsFromGeostrophicModes(kMode, lMode, jMode, phi, u);
self.A0(kIndex(abs(A0Amp)>0),lIndex(abs(A0Amp)>0),jIndex(abs(A0Amp)>0)) = A0Amp(abs(A0Amp)>0) + self.A0(kIndex(abs(A0Amp)>0),lIndex(abs(A0Amp)>0),jIndex(abs(A0Amp)>0));

self.A0 = WVTransform.makeHermitian(self.A0);

% When we hand back the actual frequency and wavenumbers, but we honor the
% users original intent and the match the signs they provided.
kMode(kMode<0) = kMode(kMode<0) + self.Nx;
lMode(lMode<0) = lMode(lMode<0) + self.Ny;
linearIndices = sub2ind(size(self.A0),kMode+1,lMode+1,jIndex);
[K,L,~] = self.kljGrid;
k = K(linearIndices);
l = L(linearIndices);
end