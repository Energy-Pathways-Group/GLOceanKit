function [omega,k,l] = SetGriddedWavesWithWavemodes(self, kMode, lMode, jMode, phi, Amp, signs)
self.RemoveAllGriddedWaves();
[omega,k,l] = self.AddGriddedWavesWithWavemodes(kMode, lMode, jMode, phi, Amp, signs);
end