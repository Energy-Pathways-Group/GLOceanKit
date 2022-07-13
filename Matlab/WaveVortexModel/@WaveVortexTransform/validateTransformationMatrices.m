function [C11,C21,C31,C12,C22,C32,C13,C23,C33] = validateTransformationMatrices(self)
% used to confirm if $$S$$ and $$S^{-1}$$ are inverses
% 
% - Topic: Validation and internal unit testing
% - Declaration: [C11,C21,C31,C12,C22,C32,C13,C23,C33] = validateTransformationMatrices(self)
% 
% This is S*S^{-1} and therefore returns the values in
% wave-vortex space. So, C11 represents Ap and should be 1s
% where we expected Ap solutions to exist.
%
% Maybe check that max(abs(C12(:))) is very small.

C11 = self.ApU.*self.UAp + self.ApV.*self.VAp + self.ApN.*self.NAp;
C21 = self.AmU.*self.UAp + self.AmV.*self.VAp + self.AmN.*self.NAp;
C31 = self.A0U.*self.UAp + self.A0V.*self.VAp + self.A0N.*self.NAp;

C12 = self.ApU.*self.UAm + self.ApV.*self.VAm + self.ApN.*self.NAm;
C22 = self.AmU.*self.UAm + self.AmV.*self.VAm + self.AmN.*self.NAm;
C32 = self.A0U.*self.UAm + self.A0V.*self.VAm + self.A0N.*self.NAm;

C13 = self.ApU.*self.UA0 + self.ApV.*self.VA0 + self.ApN.*self.NA0;
C23 = self.AmU.*self.UA0 + self.AmV.*self.VA0 + self.AmN.*self.NA0;
C33 = self.A0U.*self.UA0 + self.A0V.*self.VA0 + self.A0N.*self.NA0;
end