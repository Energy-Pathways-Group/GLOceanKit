function S_f = spectrumWithGgTransform(self,f)
arguments
    self WVDiagnostics
    f
end
wvt = self.wvt;

prefactorK = 2*ones(1,wvt.Nkl); prefactorK(1) = 1;
prefactor = wvt.g * prefactorK;

f_bar = wvt.transformFromSpatialDomainWithGg(wvt.transformFromSpatialDomainWithFourier(f));
S_f = prefactor .* abs(f_bar).^2;
end