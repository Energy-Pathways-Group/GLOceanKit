function dims = requiredDimensions(self)
% retrieve a list of required dimensions
%
% - Topic: Utility function — Metadata
arguments
    self PMAnnotatedClass {mustBeNonempty}
end
className = class(self);
dims = feval(strcat(className,'.classRequiredDimensions'));
end