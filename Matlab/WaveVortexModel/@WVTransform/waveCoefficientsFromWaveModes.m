function [kIndex,lIndex,jIndex,ApAmp,AmAmp] = waveCoefficientsFromWaveModes(self, kMode, lMode, jMode, phi, u, signs)
% Returns the indices (and re-normalized values) of the wave mode appropriate for the Ap, Am matrices.
%
% Returns the indices (and re-normalized values) of the wave mode
% appropriate for the Ap, Am matrices. This works in conjunction with the
% makeHermitian function, which then sets the appropriate conjugate. At the
% moment we made the (perhaps bad) choice that the negative l components
% are redundant, but to take advantage of the FFT, we may change this in
% the future.
% 
% For example, wave mode with l<0, is equivalent to a wave mode with l>0
% and the signs flipped on all the other quantities.
%
% The values given must meet the following requirements:
% (k0 > -Nx/2 && k0 < Nx/2)
% (l0 > -Ny/2 && l0 < Ny/2)
% (j0 >= 1 && j0 <= nModes)
% phi is in radians, from 0-2pi
% Amp is the fluid velocity U
% sign is +/-1, indicating the sign of the frequency.
%
% - Topic: Initial conditions — Waves
% - Declaration: [kIndex,lIndex,jIndex,ApAmp,AmAmp] = waveCoefficientsFromWaveModes(kMode, lMode, jMode, phi, u, signs)
% - Parameter kMode: integer index, (k0 > -Nx/2 && k0 < Nx/2)
% - Parameter lMode: integer index, (l0 > -Ny/2 && l0 < Ny/2)
% - Parameter jMode: integer index, (j0 >= 1 && j0 <= nModes), unless k=l=0 in which case j=0 is okay (inertial oscillations)
% - Parameter phi: phase in radians, (0 <= phi <= 2*pi)
% - Parameter Amp: fluid velocity u (m/s)
% - Parameter sign: sign of the frequency, +1 or -1
arguments
    self WVTransform {mustBeNonempty}
    kMode (:,1) double
    lMode (:,1) double
    jMode (:,1) double
    phi (:,1) double
    u (:,1) double
    signs (:,1) double
end

if ~isequal(size(kMode), size(lMode), size(jMode), size(phi), size(u), size(signs))
    error('All input array must be of equal size');
end

kNorm = kMode;
lNorm = lMode;
jNorm = jMode;
phiNorm = phi;
uNorm = u;
signIndex = signs;
ApAmp = zeros(size(kMode));
AmAmp = zeros(size(kMode));
for iMode = 1:length(kMode)
    % User input sanity checks. We don't deal with the Nyquist.
    if (kNorm(iMode) <= -self.Nx/2 || kNorm(iMode) >= self.Nx/2)
        error('Invalid choice for k0 (%d). Must be an integer %d < k0 < %d',kNorm(iMode),-self.Nx/2+1,self.Nx/2-1);
    end
    if (lNorm(iMode) <= -self.Ny/2 || lNorm(iMode) >= self.Ny/2)
        error('Invalid choice for l0 (%d). Must be an integer %d < l0 < %d',lNorm(iMode),-self.Ny/2+1,self.Ny/2+1);
    end

    if ~( (jNorm(iMode) == 0 && lNorm(iMode) == 0 && kNorm(iMode) == 0) || (jNorm(iMode) >= 1 && jNorm(iMode) <= self.Nz-1) )
        warning('Invalid choice for j0 (%d). Must be an integer 1 <= j <= %d, unless k=l=0, in which case j=0 is okay.',jNorm(iMode), self.Nz-1);
    end

    % Deal with the negative wavenumber cases (and inertial)
    if lNorm(iMode) == 0 && kNorm(iMode) == 0 % inertial
        if signIndex(iMode) < 1
            signIndex(iMode) = 1;
            phiNorm(iMode)  = -phiNorm;
        end
    elseif lNorm(iMode) == 0 && kNorm(iMode) < 0
        kNorm(iMode) = -kNorm(iMode);
        signIndex(iMode) = -1*signIndex(iMode);
        uNorm(iMode) = -1*uNorm(iMode);
        phiNorm(iMode) = -phiNorm(iMode);
    elseif lNorm(iMode) < 0
        lNorm(iMode) = -lNorm(iMode);
        kNorm(iMode) = -kNorm(iMode);
        signIndex(iMode) = -1*signIndex(iMode);
        uNorm(iMode) = -1*uNorm(iMode);
        phiNorm(iMode) = -phiNorm(iMode);
    end

    % Rewrap (k0,l0) to follow standard FFT wrapping. l0 should
    % already be correct.
    if (kNorm(iMode) < 0)
        kNorm(iMode) = self.Nx + kNorm(iMode);
    end

    ratio = self.uMaxGNormRatioForWave(kNorm(iMode), lNorm(iMode), jNorm(iMode));
    A = uNorm*ratio/2*exp(sqrt(-1)*phiNorm);
    if lNorm(iMode) == 0 && kNorm(iMode) == 0
        ApAmp(iMode) = A;
        AmAmp(iMode) = conj(A);
    elseif signIndex(iMode) > 0
        ApAmp(iMode) = A;
        AmAmp(iMode) = 0;
    else
        ApAmp(iMode) = 0;
        AmAmp(iMode) = A;
    end

end

kIndex = kNorm+1;
lIndex = lNorm+1;
jIndex = jNorm+1;

end