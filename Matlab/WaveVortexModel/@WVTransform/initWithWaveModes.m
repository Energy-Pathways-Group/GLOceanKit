function [omega,k,l] = initWithWaveModes(self, kMode, lMode, jMode, phi, Amp, signs)
% initialize with the given wave modes
%
% Clears variables Ap,Am,A0 and then sets the given wave modes.
% - Topic: Initial conditions — Waves
% - Declaration: [omega,k,l] = initWithWaveModes( kMode, lMode, jMode, phi, Amp, signs)
% - Parameter kMode: integer index, (k0 > -Nx/2 && k0 < Nx/2)
% - Parameter lMode: integer index, (l0 > -Ny/2 && l0 < Ny/2)
% - Parameter jMode: integer index, (j0 >= 1 && j0 <= nModes), unless k=l=0 in which case j=0 is okay (inertial oscillations)
% - Parameter phi: phase in radians, (0 <= phi <= 2*pi)
% - Parameter Amp: fluid velocity u (m/s)
% - Parameter sign: sign of the frequency, +1 or -1
% - Returns omega: frequencies of the waves (radians/s)
% - Returns k: wavenumber k of the waves (radians/m)
% - Returns l: wavenumber l of the waves (radians/m)
self.Ap = zeros(size(self.Ap));
self.Am = zeros(size(self.Am));
self.A0 = zeros(size(self.A0));
[omega,k,l] = self.setWaveModes(kMode, lMode, jMode, phi, Amp, signs);
end