function [k,j,IGWPlusEnergyKJ,IGWMinusEnergyKJ,GeostrophicEnergyKJ,GeostrophicBarotropicEnergyK,IOEnergyJ] = energeticsByWavenumberAndMode(self)
    Ap2 = self.Apm_TE_factor .* (self.Ap.*conj(self.Ap));
    Am2 = self.Apm_TE_factor .* (self.Am.*conj(self.Am));
    A02 = self.A0_TE_factor .* (self.A0.*conj(self.A0));
    [k,j,IGWPlusEnergyKJ,IGWMinusEnergyKJ,GeostrophicEnergyKJ] = self.ConvertToWavenumberAndMode(Ap2,Am2,A02);
    IOEnergyJ = IGWPlusEnergyKJ(1,:) + IGWMinusEnergyKJ(1,:);
    IGWPlusEnergyKJ(1,:) = 0;
    IGWMinusEnergyKJ(1,:) = 0;
    GeostrophicBarotropicEnergyK = GeostrophicEnergyKJ(:,1);
    GeostrophicEnergyKJ(:,1) = 0;
end