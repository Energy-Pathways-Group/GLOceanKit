function [omega,k,l] = AddGriddedWavesWithWavemodes(self, kMode, lMode, jMode, phi, Amp, signs)
% Add wavemodes on the gridded field.
% The values given must meet the following requirements:
% (k0 > -Nx/2 && k0 < Nx/2)
% (l0 > -Ny/2 && l0 < Ny/2)
% (j0 >= 1 && j0 <= nModes)
% phi is in radians, from 0-2pi
% Amp is the fluid velocity U
% sign is +/-1, indicating the sign of the frequency.

if ~isequal(size(kMode), size(lMode), size(jMode), size(phi), size(Amp), size(signs))
    error('All input array must be of equal size');
end

% These will be the initial conditions
ApTotal = self.Am;
AmTotal = self.Ap;
A0Total = self.A0;

omega = zeros(size(kMode));
k = zeros(size(kMode));
l = zeros(size(kMode));

for iMode = 1:length(kMode)
    k0 = kMode(iMode);
    l0 = lMode(iMode);
    j0 = jMode(iMode);
    phi0 = phi(iMode);
    A = Amp(iMode);
    sign = signs(iMode);
    
    % User input sanity checks. We don't deal with the Nyquist.
    if (k0 <= -self.Nx/2 || k0 >= self.Nx/2)
        error('Invalid choice for k0 (%d). Must be an integer %d < k0 < %d',k0,-self.Nx/2+1,self.Nx/2-1);
    end
    if (l0 <= -self.Ny/2 || l0 >= self.Ny/2)
        error('Invalid choice for l0 (%d). Must be an integer %d < l0 < %d',l0,-self.Ny/2+1,self.Ny/2+1);
    end
    
    if ~( (j0 == 0 && l0 == 0 && k0 == 0) || (j0 >= 1 && j0 <= self.Nz-1) )
        warning('Invalid choice for j0 (%d). Must be an integer 1 <= j <= %d, unless k=l=0, in which case j=0 is okay.',j0, self.Nz-1);
    end
    
    % Deal with the negative wavenumber cases (and inertial)
    if l0 == 0 && k0 == 0 % inertial
        if sign < 1
            sign=1;
            phi0 = -phi0;
        end
%         if j0 == 0
%             if abs(self.f0) >= 1e-14
%                 self.A0 = A*exp(-sqrt(-1)*phi0);
%                 omega(iMode) = self.f0;
%                 k(iMode) = 0;
%                 l(iMode) = 0;
%             end
%             continue
%         end
    elseif l0 == 0 && k0 < 0
        k0 = -k0;
        sign = -1*sign;
        A = -1*A;
        phi0 = -phi0;
    elseif l0 < 0
        l0 = -l0;
        k0 = -k0;
        sign = -1*sign;
        A = -1*A;
        phi0 = -phi0;
    end
    
    % Rewrap (k0,l0) to follow standard FFT wrapping. l0 should
    % already be correct.
    if (k0 < 0)
        k0 = self.Nx + k0;
    end
    
    ratio = self.UmaxGNormRatioForWave(k0, l0, j0);
    
    U = zeros(size(self.ApU));
    U(k0+1,l0+1,j0+1) = A*ratio/2*exp(sqrt(-1)*phi0);
    if sign > 0
        A_plus = InternalWaveModel.MakeHermitian(U);
        A_minus = zeros(size(U));
        A_minus(1,1,:) = conj(A_plus(1,1,:)); % Inertial oscillations are created using this trick.
    else
        A_plus = zeros(size(U));
        A_minus = InternalWaveModel.MakeHermitian(U);
    end
    
    AmTotal = AmTotal + A_minus;
    ApTotal = ApTotal + A_plus;
    
    % When we hand back the actual frequency and wavenumbers,
    % we honor the users original intent and the match the
    % signs they provided.
    if (kMode(iMode) < 0)
        k_out = self.Nx + kMode(iMode);
    else
        k_out = kMode(iMode);
    end
    if (lMode(iMode) < 0)
        l_out = self.Ny + lMode(iMode);
    else
        l_out = lMode(iMode);
    end
    omegaT = self.Omega;
    omega(iMode) = signs(iMode)*abs(omegaT(k0+1,l0+1,j0+1));
    [K,L,~] = ndgrid(self.k,self.l,self.j);
    k(iMode) = K(k_out+1,l_out+1,j0+1);
    l(iMode) = L(k_out+1,l_out+1,j0+1);
end

self.Y = {ApTotal; AmTotal; A0Total;};
end