function vars = requiredProperties(self)
% retrieve a list of required dimensions
%
% - Topic: Utility function — Metadata
arguments
    self PMAnnotatedClass {mustBeNonempty}
end
className = class(self);
vars = feval(strcat(className,'.classRequiredProperties'));
end