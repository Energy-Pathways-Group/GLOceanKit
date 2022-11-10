function [C11,C21,C31,C12,C22,C32,C13,C23,C33] = validateTransformationMatrices(self)
% used to confirm if $$S$$ and $$S^{-1}$$ are inverses
% 
% - Topic: Validation and internal unit testing
% - Declaration: [C11,C21,C31,C12,C22,C32,C13,C23,C33] = validateTransformationMatrices(self)
% 
% This is S^{-1}*S and therefore returns the values in
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

[IO,SGW,IGW,MDA,SG,IG] = self.masksForAllFlowConstituents();

fprintf('- This is a test to confirm that S^{-1} * S returns the identity matrix. This is essentially a 3x3 matrix for each (k,l,j).\n')
fprintf('- Each group of solutions is tested to confirm that it returns one along the diagonal and zero on the off-diagonal.\n')
fprintf('- Note that if a solution fails to return 1 along the diagonal, it is likely the transformation does not support that solution.\n\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inertial oscillations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yesMask = (C11 - ones(size(C11))) .* IO;
noMask = C21 .* IO + C31 .* IO;
yesError = max(abs(yesMask(:)));
noError = max(abs(noMask(:)));
if (yesError < 1e-15 && noError < 1e-15)
    fprintf('transformation matrices invert for inertial oscillations in Ap\n');
else
    fprintf('transformation matrices FAIL to invert for inertial oscillations Ap\n');
    fprintf('\t->diagonal is 1 to 1 part in 10^%d\n',round(log10(yesError)));
    fprintf('\t->off-diagonal is 0 to 1 part in 10^%d\n',round(log10(noError)));
end

yesMask = (C22 - ones(size(C22))) .* IO;
noMask = C12 .* IO + C32 .* IO;
yesError = max(abs(yesMask(:)));
noError = max(abs(noMask(:)));
if (yesError < 1e-15 && noError < 1e-15)
    fprintf('transformation matrices invert for inertial oscillations in Am\n');
else
    fprintf('transformation matrices FAIL to invert for inertial oscillationsAm\n');
    fprintf('\t->diagonal is 1 to 1 part in 10^%d\n',round(log10(yesError)));
    fprintf('\t->off-diagonal is 0 to 1 part in 10^%d\n',round(log10(noError)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% surface gravity waves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yesMask = (C11 - ones(size(C11))) .* SGW;
noMask = C21 .* SGW + C31 .* SGW;
yesError = max(abs(yesMask(:)));
noError = max(abs(noMask(:)));
if (yesError < 1e-15 && noError < 1e-15)
    fprintf('transformation matrices invert for surface gravity waves in Ap\n');
else
    fprintf('transformation matrices FAIL to invert for surface gravity waves Ap\n');
    fprintf('\t->diagonal is 1 to 1 part in 10^%d\n',round(log10(yesError)));
    fprintf('\t->off-diagonal is 0 to 1 part in 10^%d\n',round(log10(noError)));
end

yesMask = (C22 - ones(size(C22))) .* SGW;
noMask = C12 .* SGW + C32 .* SGW;
yesError = max(abs(yesMask(:)));
noError = max(abs(noMask(:)));
if (yesError < 1e-15 && noError < 1e-15)
    fprintf('transformation matrices invert for surface gravity waves in Am\n');
else
    fprintf('transformation matrices FAIL to invert for surface gravity waves Am\n');
    fprintf('\t->diagonal is 1 to 1 part in 10^%d\n',round(log10(yesError)));
    fprintf('\t->off-diagonal is 0 to 1 part in 10^%d\n',round(log10(noError)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inertial gravity waves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yesMask = (C11 - ones(size(C11))) .* IGW;
noMask = C21 .* IGW + C31 .* IGW;
yesError = max(abs(yesMask(:)));
noError = max(abs(noMask(:)));
if (yesError < 1e-15 && noError < 1e-15)
    fprintf('transformation matrices invert for internal gravity waves in Ap\n');
else
    fprintf('transformation matrices FAIL to invert for internal gravity wavesAp\n');
    fprintf('\t->diagonal is 1 to 1 part in 10^%d\n',round(log10(yesError)));
    fprintf('\t->off-diagonal is 0 to 1 part in 10^%d\n',round(log10(noError)));
end

yesMask = (C22 - ones(size(C22))) .* IGW;
noMask = C12 .* IGW + C32 .* IGW;
yesError = max(abs(yesMask(:)));
noError = max(abs(noMask(:)));
if (yesError < 1e-15 && noError < 1e-15)
    fprintf('transformation matrices invert for internal gravity waves in Am\n');
else
    fprintf('transformation matrices FAIL to invert for internal gravity wavesAm\n');
    fprintf('\t->diagonal is 1 to 1 part in 10^%d\n',round(log10(yesError)));
    fprintf('\t->off-diagonal is 0 to 1 part in 10^%d\n',round(log10(noError)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mean density anomaly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yesMask = (C33 - ones(size(C33))) .* MDA;
noMask = (C13+C23) .* MDA;
yesError = max(abs(yesMask(:)));
noError = max(abs(noMask(:)));
if (yesError < 1e-15 && noError < 1e-15)
    fprintf('transformation matrices invert for mean density anomaly solutions\n');
else
    fprintf('transformation matrices FAIL to invert for mean density anomaly solutions\n');
    fprintf('\t->diagonal is 1 to 1 part in 10^%d\n',round(log10(yesError)));
    fprintf('\t->off-diagonal is 0 to 1 part in 10^%d\n',round(log10(noError)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% surface geostrophic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yesMask = (C33 - ones(size(C33))) .* SG;
noMask = (C13+C23) .* SG;
yesError = max(abs(yesMask(:)));
noError = max(abs(noMask(:)));
if (yesError < 1e-15 && noError < 1e-15)
    fprintf('transformation matrices invert for surface geostrophic solutions\n');
else
    fprintf('transformation matrices FAIL to invert for surface geostrophic solutions\n');
    fprintf('\t->diagonal is 1 to 1 part in 10^%d\n',round(log10(yesError)));
    fprintf('\t->off-diagonal is 0 to 1 part in 10^%d\n',round(log10(noError)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interior geostrophic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yesMask = (C33 - ones(size(C33))) .* IG;
noMask = (C13+C23) .* IG;
yesError = max(abs(yesMask(:)));
noError = max(abs(noMask(:)));
if (yesError < 1e-15 && noError < 1e-15)
    fprintf('transformation matrices invert for internal geostrophic solutions\n');
else
    fprintf('transformation matrices FAIL to invert for internal geostrophic solutions\n');
    fprintf('\t->diagonal is 1 to 1 part in 10^%d\n',round(log10(yesError)));
    fprintf('\t->off-diagonal is 0 to 1 part in 10^%d\n',round(log10(noError)));
end


end