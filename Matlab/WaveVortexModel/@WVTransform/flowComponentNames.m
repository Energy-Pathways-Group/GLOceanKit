function names = flowComponentNames(self)
% retrieve the names of all available variables
%
% - Topic: Utility function — Metadata
arguments (Input)
    self WVTransform {mustBeNonempty}
end
arguments (Output)
    names string
end
names = self.flowComponentNameMap.keys;
end