% k is the frequency, is radians per meter
% d is the distance from center of mass, in meters
% z is the depth, in meters
% rho is the density
% f0 is the coriolis frequency to be used
% quick_and_dirty is a flag indicating whether or not you want be really
% careful about getting the wavelengths exactly correct. Basically, yes,
% you should use the quick and dirty method and save yourself some time.
function [Gamma] = HKEVerticalStructureFunction( z, rho, latitude )

% lets only use half the modes
nModes = floor(length(z)/2);

% For each frequency, we now try to determine the associated wavelengths
% and eigendepths
% Note that if you switch this to 'rigid_lid', then the counting of the
% modes will change. The mode in the first position will be the first
% baroclinic mode, rather than the barotropic mode.

[F, ~, h, ~] = InternalWaveModesFromDensityProfile_Spectral( rho, z, z, 0, latitude, 'total_energy', 'rigid_lid', 'fixed_k' );
h(find(h==0.0))=1e-30;

% Create the modal energy fall-off function
j_star = 3;
j=1:nModes; % [1 nModes]
H = (j_star)^(1.5)./(j+j_star).^(5/2);
%normalization = sum(H);
normalization = 2/3;
H = H./normalization;

% part of the structure function definition
H = H./h(1:nModes);

% size(F2)=[length(z) nModes]
F2 = F(:,1:nModes).*F(:,1:nModes);

Gamma = sum(repmat(H,[length(z) 1]).*F2,2);