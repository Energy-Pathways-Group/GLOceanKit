function [omega,k,l] = addWaveModes(self, kMode, lMode, jMode, phi, u, signs)
[kIndex,lIndex,jIndex,ApAmp,AmAmp] = self.waveCoefficientsFromWaveModes(kMode, lMode, jMode, phi, u, signs);
self.Ap(kIndex(abs(ApAmp)>0),lIndex(abs(ApAmp)>0),jIndex(abs(ApAmp)>0)) = ApAmp(abs(ApAmp)>0) + self.Ap(kIndex(abs(ApAmp)>0),lIndex(abs(ApAmp)>0),jIndex(abs(ApAmp)>0));
self.Am(kIndex(abs(AmAmp)>0),lIndex(abs(AmAmp)>0),jIndex(abs(AmAmp)>0)) = AmAmp(abs(AmAmp)>0) + self.Am(kIndex(abs(AmAmp)>0),lIndex(abs(AmAmp)>0),jIndex(abs(AmAmp)>0));

self.Ap = WaveVortexTransform.MakeHermitian(self.Ap);
self.Am = WaveVortexTransform.MakeHermitian(self.Am);

% When we hand back the actual frequency and wavenumbers, but we honor the
% users original intent and the match the signs they provided.
kMode(kMode<0) = kMode(kMode<0) + self.Nx;
lMode(lMode<0) = lMode(lMode<0) + self.Ny;
linearIndices = sub2ind(size(Matrix),kMode+1,lMode+1,jIndex);
omegaT = self.Omega;
omega = signs.*abs(omegaT(linearIndices));
[K,L,~] = ndgrid(self.k,self.l,self.j);
k = K(linearIndices);
l = L(linearIndices);
end