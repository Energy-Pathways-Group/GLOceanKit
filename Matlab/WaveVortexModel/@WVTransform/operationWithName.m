function val = operationWithName(self,name)
% retrieve a WVOperation by name
%
% - Topic: Utility function — Metadata
arguments
    self WVTransform {mustBeNonempty}
    name char {mustBeNonempty}
end
val = self.operationNameMap(name);
end