%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Returns a sets of 'masks' (matrices with 1s or 0s) indicating where
% different solution types live in the Ap, Am, A0 matrices.
%
% IO and IGW indicate where the inertial oscillation (IO) and internal
% gravity wave (IGW) solutions live in the Ap and Am matrices.
%
% G, G0, and Rhobar indicate where the geostrophic (G), barotropic
% geostrophic (G0) and mean density anomaly (Rhobar) solutions live in the
% A0 matrix.
% 
% For example, if you define A = IGW .* Ap; then A will contain only the
% positive frequency internal gravity solutions.
function [IO,IGW,G,G0,Rhobar] = FlowConstituentMasks(self)
    % inertial oscillations only exist at k=l=0
    IO = zeros(size(self.Ap));
    IO(1,1,:) = 1;

    % IGWs, zero out all j=0 (no solutions), and k=l=0 values (inertial).
    IGW = ones(size(self.Ap));
    IGW(:,:,1) = 0;
    IGW(1,1,:) = 0;

    % barotropic geostrophic exists at all k,l>0, j=0.
    G0 = zeros(size(self.Ap));
    G0(:,:,1) = 1;
    G0(1,1,1) = 0;

    % mean density anomaly exists at k=l=0, j>1.
    Rhobar = zeros(size(self.Ap));
    Rhobar(1,1,2:end) = 1;

    % zero out all j=0, and k=l=0 values.
    G = ones(size(self.Ap));
    G(1,1,:) = 0;
    G(:,:,1) = 0;
end