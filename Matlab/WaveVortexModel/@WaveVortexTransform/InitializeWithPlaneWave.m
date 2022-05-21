function period = InitializeWithPlaneWave(self, k0, l0, j0, UAmp, sign)
omega = self.SetGriddedWavesWithWavemodes(k0,l0,j0,0,UAmp,sign);
period = 2*pi/abs(omega);
end