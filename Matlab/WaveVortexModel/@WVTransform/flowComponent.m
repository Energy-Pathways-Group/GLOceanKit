function val = flowComponent(self,name)
% retrieve a WVFlowComponent by name
%
% - Topic: Utility function — Metadata
arguments
    self WVTransform {mustBeNonempty}
    name char {mustBeNonempty}
end
val = self.flowComponentNameMap(name);
end