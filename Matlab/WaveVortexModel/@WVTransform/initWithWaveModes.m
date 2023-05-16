function [omega,k,l] = initWithWaveModes(self, waveproperties)
% initialize with the given wave modes
%
% $$
% sin(k*x+l*y)*F_j*sin(omega*t + phi)
% $$
%
% Clears variables Ap,Am,A0 and then sets the given wave modes.
% - Topic: Initial conditions — Waves
% - Declaration: [omega,k,l] = initWithWaveModes( kMode, lMode, jMode, phi, Amp, signs)
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
    waveproperties.k (1,1) double
    waveproperties.l (1,1) double
    waveproperties.j (1,1) double
    waveproperties.phi (1,1) double
    waveproperties.u (1,1) double
    waveproperties.sign (1,1) double
end
self.Ap = zeros(size(self.Ap));
self.Am = zeros(size(self.Am));
self.A0 = zeros(size(self.A0));

[omega,k,l] = self.setWaveModes(k=waveproperties.k,l=waveproperties.l,j=waveproperties.j,phi=waveproperties.phi,u=waveproperties.u,sign=waveproperties.sign);

end
