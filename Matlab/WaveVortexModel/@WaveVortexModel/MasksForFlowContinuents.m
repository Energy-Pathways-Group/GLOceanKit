function [ApmMask,A0Mask] = MasksForFlowContinuents(self,flowConstituents)

ApmMask = zeros(size(self.Ap));
A0Mask = zeros(size(self.Ap));

if bitand(flowConstituents.inertialOscillation,flowConstituents)
    IO = zeros(size(self.Ap));
    IO(1,1,:) = 1;

    ApmMask = ApmMask || IO;
end

end