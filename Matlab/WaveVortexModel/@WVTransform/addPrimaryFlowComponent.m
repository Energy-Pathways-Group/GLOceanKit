function addPrimaryFlowComponent(self,primaryFlowComponent)
% add a primary flow component, automatically added to the flow
% components
%
% - Topic: Utility function — Metadata
arguments
    self WVTransform {mustBeNonempty}
    primaryFlowComponent (1,:) WVPrimaryFlowComponent {mustBeNonempty}
end
for i=1:length(primaryFlowComponent)
    self.primaryFlowComponentNameMap(primaryFlowComponent(i).shortName) = primaryFlowComponent(i);
    self.flowComponentNameMap(primaryFlowComponent(i).shortName) = primaryFlowComponent(i);
end
end