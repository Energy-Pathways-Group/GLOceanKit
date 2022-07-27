function [ApmMask,A0Mask] = masksForFlowConstituents(self,flowConstituents)
% Returns a sets of 'masks' indicating where different solution types live in the Ap, Am, A0 matrices.
%
% Returns a sets of 'masks' (matrices with 1s or 0s) indicating where
% different solution types live in the Ap, Am, A0 matrices.
%
% Basic usage,
% ```matlab
% [ApmMask,A0Mask] = wvm.masksForFlowConstituents(WVFlowConstituents('internalGravityWave','inertialOscillation');
% ```
% will return a mask that contains 1 at the locations of the internal
% gravity waves and inertial oscillations in the Ap/Am matrices. All other
% entries will be zero.
%
% For example, if you define ``A = ApmMask .* Ap;`` then A will contain only the
% positive frequency internal gravity solutions and half the inertial solutions.
%
% - Topic: Masks
% - Declaration: [ApmMask,A0Mask] = masksForFlowConstituents(flowConstituents)
% - Parameter flowConstituents: `WVFlowConstituents` type
% - Returns ApmMask: mask for the Ap and Am matrices
% - Returns A0Mask: mask for the A0 matrix
ApmMask = zeros(size(self.Ap));
A0Mask = zeros(size(self.A0));

[IO,SGW,IGW,MDA,SG,IG] = self.masksForAllFlowConstituents();

if flowConstituents.Contains(WVFlowConstituents.inertialOscillation)
    ApmMask = ApmMask | IO;
end

if flowConstituents.Contains(WVFlowConstituents.surfaceGravityWave)
    ApmMask = ApmMask | SGW;
end

if flowConstituents.Contains(WVFlowConstituents.internalGravityWave)
    ApmMask = ApmMask | IGW;
end

if flowConstituents.Contains(WVFlowConstituents.meanDensityAnomaly)
    A0Mask = A0Mask | MDA;
end

if flowConstituents.Contains(WVFlowConstituents.surfaceGeostrophic)
    A0Mask = A0Mask | SG;
end

if flowConstituents.Contains(WVFlowConstituents.internalGeostrophic)
    A0Mask = A0Mask | IG;
end

end