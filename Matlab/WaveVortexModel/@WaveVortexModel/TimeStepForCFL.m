function [deltaT,advectiveDT,oscillatoryDT] = TimeStepForCFL(self, cfl, outputInterval)
% Return the velocity field, which is the sum of the gridded
% and external/free waves at time t. Note that if you do not
% need w, don't request it and it won't be computed.
if nargin == 1
    cfl = 0.25;
end

omega = self.Omega;
period = 2*pi/max(abs(omega(:)));
[u,v] = self.VelocityFieldAtTime(0.0);
U = max(max(max( sqrt(u.*u + v.*v) )));
dx = (self.x(2)-self.x(1));

advectiveDT = cfl*dx/U;
oscillatoryDT = cfl*period;
% A cfl of 1/12 for oscillatoryDT might be necessary for good numerical precision when advecting particles.

fprintf('dX/U = %.1f s (%.1f min). The highest frequency resolved IGW has period of %.1f s (%.1f min).\n', dx/U,dx/U/60,period,period/60);

if advectiveDT < oscillatoryDT
    deltaT = advectiveDT;
else
    deltaT = oscillatoryDT;
end

if nargin == 3
    deltaT = outputInterval/ceil(outputInterval/deltaT);
    stepsPerOutput = round(outputInterval/deltaT);
    fprintf('Rounding to match the output interval dt: %.2f s (%d steps per output)\n',deltaT,stepsPerOutput);
end

end