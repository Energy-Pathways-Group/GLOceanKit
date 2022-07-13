function removeAllInertialMotions(self)
% remove all inertial motions
%
% All inertial motions are removed
% - Topic: Initial conditions — Inertial Oscillations
% - Declaration: removeAllInertialMotions()
self.Ap(1,1,:) = 0;
self.Am(1,1,:) = 0;
end