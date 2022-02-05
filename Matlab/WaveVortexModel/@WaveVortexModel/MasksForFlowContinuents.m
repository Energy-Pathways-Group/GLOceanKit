%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Returns a sets of 'masks' (matrices with 1s or 0s) indicating where
% different solution types live in the Ap, Am, A0 matrices.
%
% Basic usage,
% [ApmMask,A0Mask] = wvm.MasksForFlowContinuents(FlowConstituents('internalGravityWave','inertialOscillation');
% will return a mask that contains 1 at the locations of the internal
% gravity waves and inertial oscillations in the Ap/Am matrices. All other
% entries will be zero.
%
% For example, if you define A = ApmMask .* Ap; then A will contain only the
% positive frequency internal gravity solutions and half the inertial solutions.
function [ApmMask,A0Mask] = MasksForFlowContinuents(self,flowConstituents)

ApmMask = zeros(size(self.Ap));
A0Mask = zeros(size(self.A0));

[IO,SGW,IGW,MDA,SG,IG] = self.MasksForAllFlowConstituents();

if flowConstituents.Contains(FlowConstituents.inertialOscillation)
    ApmMask = ApmMask | IO;
elseif flowConstituents.Contains(FlowConstituents.surfaceGravityWave)
    ApmMask = ApmMask | SGW;
elseif flowConstituents.Contains(FlowConstituents.internalGravityWave)
    ApmMask = ApmMask | IGW;
elseif flowConstituents.Contains(FlowConstituents.meanDensityAnomaly)
    A0Mask = A0Mask | MDA;
elseif flowConstituents.Contains(FlowConstituents.surfaceGeostrophic)
    A0Mask = A0Mask | SG;
elseif flowConstituents.Contains(FlowConstituents.internalGeostrophic)
    A0Mask = A0Mask | IG;
end

end