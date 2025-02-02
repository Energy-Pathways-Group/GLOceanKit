function components = flowComponents(self)
arguments (Input)
    self WVTransform
end
arguments (Output)
    components WVFlowComponent
end
components = [self.flowComponentNameMap{self.flowComponentNameMap.keys}];
end