function addFlowComponent(self,flowComponent)
% add a flow component
%
% - Topic: Utility function — Metadata
arguments
    self WVTransform {mustBeNonempty}
    flowComponent (1,:) WVPrimaryFlowComponent {mustBeNonempty}
end
for i=1:length(flowComponent)
    self.flowComponentNameMap(flowComponent(i).name) = flowComponent(i);
end
end