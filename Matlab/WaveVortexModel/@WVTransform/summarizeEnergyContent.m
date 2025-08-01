function summarizeEnergyContent(self)
% displays a summary of the energy content of the fluid
%
% - Topic: Energetics
total = self.totalEnergy;
pctString = sprintf('(');
nameString = sprintf('(');
primaryFlowComponentNames_ = self.primaryFlowComponentNames;
for iVar = 1:length(primaryFlowComponentNames_)
    name = primaryFlowComponentNames_(iVar);
    pctString = [pctString,sprintf('%.1f,',100*self.totalEnergyOfFlowComponent(self.flowComponentWithName(name{1}))/total)];
    nameString = [nameString,sprintf('%s, ',name{1})];
end

fprintf('%.2g m^3/s^2 total depth integrated energy, split %s) between %s)\n',total,pctString,nameString);
end