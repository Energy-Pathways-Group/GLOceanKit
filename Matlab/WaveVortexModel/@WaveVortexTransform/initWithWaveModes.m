function [omega,k,l] = initWithWaveModes(self, kMode, lMode, jMode, phi, Amp, signs)
self.Ap = zeros(size(self.Ap));
self.Am = zeros(size(self.Am));
self.A0 = zeros(size(self.A0));
[omega,k,l] = self.setWaveModes(kMode, lMode, jMode, phi, Amp, signs);
end