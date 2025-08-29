function S_f = crossSpectrumWithGgTransform(self,phi,gamma)
arguments
    self WVGeometryDoublyPeriodicStratified
    phi
    gamma
end

prefactorK = 2*ones(1,self.Nkl); prefactorK(1) = 1;
prefactor = self.g * prefactorK;

phi_bar = self.transformFromSpatialDomainWithGg(self.transformFromSpatialDomainWithFourier(phi));
gamma_bar = self.transformFromSpatialDomainWithGg(self.transformFromSpatialDomainWithFourier(gamma));
S_f = prefactor .* real(phi_bar .* conj(gamma_bar));
end