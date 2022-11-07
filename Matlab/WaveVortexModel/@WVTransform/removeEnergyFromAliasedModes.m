function removeEnergyFromAliasedModes(self,options)
% remove all energy from aliased modes
%
% Removes the energy from modes that will alias with a quadratic
% multiplication. If running a nonlinear simulation, it is recommended that
% you perform this step before time-stepping.
% - Topic: Initial conditions
% - Declaration: removeEnergyFromAliasedModes(options)
% - Parameter jFraction: (optional) fraction of vertical mode to assume are not aliased (default 2/3)
arguments
    self WVTransform {mustBeNonempty}
    options.jFraction double {mustBePositive(options.jFraction),mustBeLessThanOrEqual(options.jFraction,1)} = 2/3
end
AntiAliasMask = self.maskForAliasedModes(jFraction=options.jFraction);
self.A0 = self.A0 .* ~AntiAliasMask;
self.Am = self.Am .* ~AntiAliasMask;
self.Ap = self.Ap .* ~AntiAliasMask;
end