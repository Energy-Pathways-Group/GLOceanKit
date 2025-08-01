function addFlowComponent(self,flowComponent)
% add a flow component
%
% - Topic: Flow components
arguments
    self WVTransform {mustBeNonempty}
    flowComponent (1,:) WVFlowComponent {mustBeNonempty}
end
for i=1:length(flowComponent)
    self.flowComponentNameMap{flowComponent(i).shortName} = flowComponent(i);
end
end