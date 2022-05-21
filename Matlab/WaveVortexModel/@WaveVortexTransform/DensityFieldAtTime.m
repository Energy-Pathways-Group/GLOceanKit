function rho = DensityFieldAtTime(self, t)
% Return the density field, which is the sum of the density
% mean field (variable in z) and the perturbation field
% (variable in time and space).
rho_bar = self.DensityMeanField;
rho_prime = self.VariableFieldsAtTime(t,'rho_prime');
rho = rho_bar + rho_prime;
end