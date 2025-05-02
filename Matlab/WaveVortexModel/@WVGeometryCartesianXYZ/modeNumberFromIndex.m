function [kMode,lMode,jMode] = modeNumberFromIndex(self,linearIndex)
arguments (Input)
    self WVTransform {mustBeNonempty}
    linearIndex (:,1) double {mustBeInteger,mustBePositive}
end
arguments (Output)
    kMode (:,1) double {mustBeInteger}
    lMode (:,1) double {mustBeInteger}
    jMode (:,1) double {mustBeInteger,mustBeNonnegative}
end
[jIndex,klIndex] = ind2sub(self.spectralMatrixSize,linearIndex);
[kMode,lMode] = self.klModeNumberFromIndex(klIndex);
jMode = jIndex - 1;
end