function summarizeEnergyContent(self)
% displays a summary of the energy content of the fluid
%
% - Topic: Energetics
total = self.totalEnergy;
pctString = sprintf('(');
nameString = sprintf('(');
for name = keys(self.primaryFlowComponentNameMap)
    pctString = [pctString,sprintf('%.1f,',100*self.totalEnergyOfFlowComponent(self.flowComponent(name{1}))/total)];
    nameString = [nameString,sprintf('%s, ',name{1})];
end

fprintf('%.2g m^3/s^2 total depth integrated energy, split %s) between %s)\n',total,pctString,nameString);
end