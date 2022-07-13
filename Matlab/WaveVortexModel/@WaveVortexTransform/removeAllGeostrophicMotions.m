function removeAllGeostrophicMotions(self)
% remove all geostrophic motions
%
% All geostrophic motions are removed by setting A0 to zero.
% - Topic: Initial conditions — Geostrophic Motions
% - Declaration: removeAllGeostrophicMotions()
self.A0 = 0*self.A0;
end