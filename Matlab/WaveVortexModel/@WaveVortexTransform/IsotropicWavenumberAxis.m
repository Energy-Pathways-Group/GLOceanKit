function kIso = IsotropicWavenumberAxis(self)
% Create a reasonable wavenumber axis
allKs = unique(reshape(abs(self.Kh),[],1),'sorted');
deltaK = max(diff(allKs));
kAxis = 0:deltaK:max(allKs);

% Thi is the final output axis for wavenumber
kIso = reshape(kAxis(1:(length(kAxis)-1)),[],1);
end