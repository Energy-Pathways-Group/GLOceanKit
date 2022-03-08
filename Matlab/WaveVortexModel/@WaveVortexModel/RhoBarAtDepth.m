function rho = RhoBarAtDepth(self,z)
    g = 9.81;
    rho = -(self.N0*self.N0*self.rho0/g)*z + self.rho0;
end
