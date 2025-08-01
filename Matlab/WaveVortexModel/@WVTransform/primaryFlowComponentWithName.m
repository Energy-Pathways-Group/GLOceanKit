function val = primaryFlowComponentWithName(self,name)
% retrieve a WVPrimaryFlowComponent by name
%
% - Topic: Flow components
arguments
    self WVTransform {mustBeNonempty}
    name char {mustBeNonempty}
end
val = self.primaryFlowComponentNameMap{name};
end