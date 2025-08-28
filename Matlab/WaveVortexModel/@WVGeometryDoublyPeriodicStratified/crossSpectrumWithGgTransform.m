function S_f = crossSpectrumWithGgTransform(self,phi,gamma)
arguments
    self WVDiagnostics
    phi
    gamma
end
wvt = self.wvt;

prefactorK = 2*ones(1,wvt.Nkl); prefactorK(1) = 1;
prefactor = wvt.g * prefactorK;

phi_bar = wvt.transformFromSpatialDomainWithGg(wvt.transformFromSpatialDomainWithFourier(phi));
gamma_bar = wvt.transformFromSpatialDomainWithGg(wvt.transformFromSpatialDomainWithFourier(gamma));
S_f = prefactor .* real(phi_bar .* conj(gamma_bar));
end