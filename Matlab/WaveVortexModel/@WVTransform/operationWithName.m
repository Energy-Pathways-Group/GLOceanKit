function val = operationWithName(self,name)
% retrieve a WVOperation by name
%
% - Topic: Utility function — Metadata
arguments (Input)
    self WVTransform {mustBeNonempty}
    name char {mustBeNonempty}
end
arguments (Output)
    val WVOperation 
end
val = self.operationNameMap{name};
end