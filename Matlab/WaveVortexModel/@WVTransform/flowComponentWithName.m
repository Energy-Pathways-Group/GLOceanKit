function val = flowComponentWithName(self,name)
% retrieve a WVFlowComponent by name
%
% - Topic: Flow components
arguments
    self WVTransform {mustBeNonempty}
    name char {mustBeNonempty}
end
val = self.flowComponentNameMap{name};
end