%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Returns six 'masks' (matrices with 1s or 0s) indicating where the six
% different solution types live in the Ap, Am, A0 matrices.
%
% IO, SGW, and IGW indicate where the inertial oscillation (IO), surface
% gravity waves (SGW), and internal gravity wave (IGW) solutions live in
% the Ap and Am matrices.
%
% MDA, SG, and IG indicate where the mean density anomaly (MDA), surface
% geostrophic (SG), and interior geostrophic (IG) solutions live in the A0
% matrix.
% 
% For example, if you define A = IGW .* Ap; then A will contain only the
% positive frequency internal gravity solutions.
function [IO,SGW,IGW,MDA,SG,IG] = MasksForAllFlowConstituents(self)
    IO = zeros(size(self.Ap));
    IO(1,1,:) = 1;

    SGW = zeros(size(self.Ap));
    SGW(:,:,1) = 1;

    IGW = ones(size(self.Ap));
    IGW(:,:,1) = 0;
    IGW(1,1,:) = 0;

    MDA = zeros(size(self.Ap));
    MDA(1,1,:) = 1;

    SG = zeros(size(self.Ap));
    SG(:,:,1) = 1;

    IG = ones(size(self.Ap));
    IG(:,:,1) = 0;
    IG(1,1,:) = 0;
end