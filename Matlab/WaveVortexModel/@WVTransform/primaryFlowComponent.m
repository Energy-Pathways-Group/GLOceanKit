function val = primaryFlowComponent(self,name)
% retrieve a WVPrimaryFlowComponent by name
%
% - Topic: Utility function — Metadata
arguments
    self WVTransform {mustBeNonempty}
    name char {mustBeNonempty}
end
val = self.primaryFlowComponentNameMap(name);
end