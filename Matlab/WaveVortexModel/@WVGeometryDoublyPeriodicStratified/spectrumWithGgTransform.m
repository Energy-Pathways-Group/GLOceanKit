function S_f = spectrumWithGgTransform(self,f)
arguments
    self WVGeometryDoublyPeriodicStratified
    f
end

prefactorK = 2*ones(1,self.Nkl); prefactorK(1) = 1;
prefactor = self.g * prefactorK;

f_bar = self.transformFromSpatialDomainWithGg(self.transformFromSpatialDomainWithFourier(f));
S_f = prefactor .* abs(f_bar).^2;
end