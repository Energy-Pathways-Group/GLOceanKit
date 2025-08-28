function S_f = spectrumWithFgTransform(self,f)
arguments
    self WVGeometryDoublyPeriodicStratified
    f
end
prefactorJ = self.h_0; prefactorJ(1) = self.Lz;
prefactorK = 2*ones(1,self.Nkl); prefactorK(1) = 1;
prefactor = prefactorJ * prefactorK;

f_bar = self.transformFromSpatialDomainWithFg(self.transformFromSpatialDomainWithFourier(f));
S_f = prefactor.*abs(f_bar).^2;
end