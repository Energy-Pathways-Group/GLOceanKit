function [omega,k,l] = setWaveModes(self, kMode, lMode, jMode, phi, u, signs)
% set amplitudes of the given wave modes
%
% Overwrite any existing wave modes with the given new values
% - Topic: Initial conditions — Waves
% - Declaration: [omega,k,l] = setWaveModes(kMode, lMode, jMode, phi, u, signs)
% - Parameter kMode: integer index, (k0 > -Nx/2 && k0 < Nx/2)
% - Parameter lMode: integer index, (l0 > -Ny/2 && l0 < Ny/2)
% - Parameter jMode: integer index, (j0 >= 1 && j0 <= nModes), unless k=l=0 in which case j=0 is okay (inertial oscillations)
% - Parameter phi: phase in radians, (0 <= phi <= 2*pi)
% - Parameter Amp: fluid velocity u (m/s)
% - Parameter sign: sign of the frequency, +1 or -1
% - Returns omega: frequencies of the waves (radians/s)
% - Returns k: wavenumber k of the waves (radians/m)
% - Returns l: wavenumber l of the waves (radians/m)
[kIndex,lIndex,jIndex,ApAmp,AmAmp] = self.waveCoefficientsFromWaveModes(kMode, lMode, jMode, phi, u, signs);
self.Ap(kIndex(abs(ApAmp)>0),lIndex(abs(ApAmp)>0),jIndex(abs(ApAmp)>0)) = ApAmp(abs(ApAmp)>0);
self.Am(kIndex(abs(AmAmp)>0),lIndex(abs(AmAmp)>0),jIndex(abs(AmAmp)>0)) = AmAmp(abs(AmAmp)>0);

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