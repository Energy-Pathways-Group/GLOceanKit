function [omega, alpha, k, l, mode, phi, A, norm] = waveModesFromWaveCoefficients(self)
% Returns normalized amplitudes and phases of all waves
%
% This returns the properties of the waves being used in the
% gridded simulation, as their properly normalized individual
% wave components. Very useful for debugging.
%
% Note that A_plus and A_minus each have half the inertial
% energy. This can be misleading, but the phasing is chosen to
% make it work. Never-the-less, we double/zero that component.
%
% - Topic: Initial conditions — Waves
% - Declaration: [omega, alpha, k, l, mode, phi, A, norm] = waveModesFromWaveCoefficients()
A_p = self.Ap;
A_p(1,1,:) = 2*A_p(1,1,:);
A_m = self.Am;
A_m(1,1,:) = 0*A_m(1,1,:);

[K,L,J] = ndgrid(self.k,self.l,self.j);
Omega = self.Omega;

[A_plus,phi_plus,linearIndex] = WaveVortexTransform.ExtractNonzeroWaveProperties(A_p);
omega_plus = Omega(linearIndex);
mode_plus = J(linearIndex);
alpha_plus = atan2(L(linearIndex),K(linearIndex));
k_plus = K(linearIndex);
l_plus = L(linearIndex);

[A_minus,phi_minus,linearIndex] = WaveVortexTransform.ExtractNonzeroWaveProperties(A_m);
omega_minus = -Omega(linearIndex);
mode_minus = J(linearIndex);
alpha_minus = atan2(L(linearIndex),K(linearIndex));
k_minus = K(linearIndex);
l_minus = L(linearIndex);

k = [k_plus; k_minus];
l = [l_plus; l_minus];
omega = [omega_plus; omega_minus];
mode = [mode_plus; mode_minus];
alpha = [alpha_plus; alpha_minus];
phi = [phi_plus; phi_minus];
A = [A_plus; A_minus];
norm = Normalization.kConstant;
end