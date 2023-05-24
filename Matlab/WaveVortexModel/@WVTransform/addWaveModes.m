function [omega,k,l] = addWaveModes(self, waveproperties)
% add amplitudes of the given wave modes
%
% Add new amplitudes to any existing amplitudes
% - Topic: Initial conditions — Waves
% - Declaration: [omega,k,l] = addWaveModes(self, kMode, lMode, jMode, phi, u, signs)
% - Parameter k: integer index, (k0 > -Nx/2 && k0 < Nx/2)
% - Parameter l: integer index, (l0 > -Ny/2 && l0 < Ny/2)
% - Parameter j: integer index, (j0 >= 1 && j0 <= nModes), unless k=l=0 in which case j=0 is okay (inertial oscillations)
% - Parameter phi: phase in radians, (0 <= phi <= 2*pi)
% - Parameter u: fluid velocity u (m/s)
% - Parameter sign: sign of the frequency, +1 or -1
% - Returns omega: frequencies of the waves (radians/s)
% - Returns k: wavenumber k of the waves (radians/m)
% - Returns l: wavenumber l of the waves (radians/m)
arguments
    self WVTransform {mustBeNonempty}
    waveproperties.k (:,1) double
    waveproperties.l (:,1) double
    waveproperties.j (:,1) double
    waveproperties.phi (:,1) double
    waveproperties.u (:,1) double
    waveproperties.sign (:,1) double
end
kMode = waveproperties.k;
lMode = waveproperties.l;
jMode = waveproperties.j;
phi = waveproperties.phi;
u = waveproperties.u;
signs = waveproperties.sign;

[kIndex,lIndex,jIndex,ApAmp,AmAmp] = self.waveCoefficientsFromWaveModes(kMode, lMode, jMode, phi, u, signs);
self.Ap(kIndex(abs(ApAmp)>0),lIndex(abs(ApAmp)>0),jIndex(abs(ApAmp)>0)) = ApAmp(abs(ApAmp)>0) + self.Ap(kIndex(abs(ApAmp)>0),lIndex(abs(ApAmp)>0),jIndex(abs(ApAmp)>0));
self.Am(kIndex(abs(AmAmp)>0),lIndex(abs(AmAmp)>0),jIndex(abs(AmAmp)>0)) = AmAmp(abs(AmAmp)>0) + self.Am(kIndex(abs(AmAmp)>0),lIndex(abs(AmAmp)>0),jIndex(abs(AmAmp)>0));

self.Ap = WVTransform.makeHermitian(self.Ap);
self.Am = WVTransform.makeHermitian(self.Am);

% When we hand back the actual frequency and wavenumbers, but we honor the
% users original intent and the match the signs they provided.
kMode(kMode<0) = kMode(kMode<0) + self.Nx;
lMode(lMode<0) = lMode(lMode<0) + self.Ny;
linearIndices = sub2ind(size(self.Ap),kMode+1,lMode+1,jIndex);
omegaT = self.Omega;
omega = signs.*abs(omegaT(linearIndices));
[K,L,~] = ndgrid(self.k,self.l,self.j);
k = K(linearIndices);
l = L(linearIndices);
end