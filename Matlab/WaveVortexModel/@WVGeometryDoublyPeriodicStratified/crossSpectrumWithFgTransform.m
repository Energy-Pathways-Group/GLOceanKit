function S_f = crossSpectrumWithFgTransform(self,phi,gamma)
arguments
    self WVGeometryDoublyPeriodicStratified
    phi
    gamma
end

prefactorJ = self.h_0; prefactorJ(1) = self.Lz;
prefactorK = 2*ones(1,self.Nkl); prefactorK(1) = 1;
prefactor = prefactorJ * prefactorK;

phi_bar = self.transformFromSpatialDomainWithFg(self.transformFromSpatialDomainWithFourier(phi));
gamma_bar = self.transformFromSpatialDomainWithFg(self.transformFromSpatialDomainWithFourier(gamma));
S_f = prefactor .* real(phi_bar .* conj(gamma_bar));
end