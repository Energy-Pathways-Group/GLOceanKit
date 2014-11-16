% the sshFD is expected to be in a wrapped, matlab format.
function [kMag, energyMag] = EnergySpectrumFromSSH( sshFD, k, l, g, f0, L_R )

% FFT algorithms wrap wavenumbers, we need to unwrap them
k = circshift( k, (length(k)/2));
l = circshift( l, (length(l)/2));
sshFD = fftshift( sshFD );

[K, L] = meshgrid(k,l);
rhoK = sqrt( K.*K + L.*L);

deltaX = 1/(2*k(length(k)/2));
deltaY = 1/(2*k(length(k)/2));
spectrumScale = deltaX*deltaY/( length(k) * length(l) );
heightPSD = spectrumScale * real( sshFD .* conj(sshFD)); % m^2 * m^2
energyPSD = ((g/f0)^2)* ( (2*2*pi*pi) * rhoK .* rhoK .* heightPSD + (1/L_R)*(1/L_R) * heightPSD ); % m^2 * m^2 / s^2

deltaK = ones(size(energyPSD)) * (l(2)-l(1));
intEnergyPSD = energyPSD .* deltaK;

% Build the wavenumber vector and find the energy in those bins
%kMag = k(length(k)/2+1:end);
%dK = kMag(2)-kMag(1);
dK = k(2)-k(1);
kMag = 0:dK:max(max(rhoK));
indices = cell( length(kMag), 1);
for i = 1:length(kMag)
	indices{i} = find( rhoK >= kMag(i)-dK/2 & rhoK < kMag(i)+dK/2 );
end

% Now sum up the energy
energyMag = zeros( size(kMag) );
for i = 1:length(kMag)
	energyMag(i) = sum( intEnergyPSD( indices{i} ) );
end