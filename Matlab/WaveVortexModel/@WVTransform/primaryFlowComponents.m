function components = primaryFlowComponents(self)
arguments (Input)
    self WVTransform
end
arguments (Output)
    components WVPrimaryFlowComponent
end
components = [self.primaryFlowComponentNameMap{self.primaryFlowComponentNameMap.keys}];
end