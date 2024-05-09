function summarizeEnergyContent(self)
% displays a summary of the energy content of the fluid
%
% - Topic: Energetics
total = self.totalEnergy;
ioPct = 100*self.inertialEnergy/total;
wavePct = 100*self.waveEnergy/total;
gPct = 100*self.geostrophicEnergy/total;
mdaPct = 100*self.mdaEnergy/total;

fprintf('%.2g m^3/s^2 total depth integrated energy, split (%.1f,%.1f,%.1f,%.1f) between (inertial,wave,geostrophic,mda)\n',total,ioPct,wavePct,gPct,mdaPct);
end