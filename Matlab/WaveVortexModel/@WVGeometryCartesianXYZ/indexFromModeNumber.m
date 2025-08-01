function index = indexFromModeNumber(self,kMode,lMode,jMode)
% return the linear index into a spectral matrix given (k,l,j)
%
% This function will return the linear index in a spectral
% matrix given a mode number.
%
% - Topic: Index Gymnastics
% - Declaration: index = indexFromModeNumber(kMode,lMode,jMode)
% - Parameter kMode: integer
% - Parameter lMode: integer
% - Returns index: a non-negative integer number
arguments (Input)
    self WVTransform {mustBeNonempty}
    kMode (:,1) double {mustBeInteger}
    lMode (:,1) double {mustBeInteger}
    jMode (:,1) double {mustBeInteger,mustBeNonnegative}
end
arguments (Output)
    index (:,1) double {mustBeInteger,mustBePositive}
end
if ~self.isValidModeNumber(kMode,lMode,jMode)
    error('Invalid WV mode number!');
end
[kMode,lMode] = self.primaryKLModeNumberFromKLModeNumber(kMode,lMode);
klIndex = self.indexFromKLModeNumber(kMode,lMode);
index = sub2ind(self.spectralMatrixSize,jMode+1,klIndex);
end