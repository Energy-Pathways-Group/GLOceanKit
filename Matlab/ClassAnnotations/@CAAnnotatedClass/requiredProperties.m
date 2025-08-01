function vars = requiredProperties(self)
% retrieve a list of required dimensions
%
% - Topic: Utility function — Metadata
arguments
    self CAAnnotatedClass {mustBeNonempty}
end
className = class(self);
vars = feval(strcat(className,'.classRequiredPropertyNames'));
end