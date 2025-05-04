function bool = isValidModeNumber(self,kMode,lMode,jMode)
% returns a boolean indicating whether (k,l,j) is a valid mode number
%
% returns a boolean indicating whether (k,l,j) is a valid mode
% number
%
% - Topic: Index Gymnastics
% - Declaration: index = isValidModeNumber(kMode,lMode,jMode)
% - Parameter kMode: integer
% - Parameter lMode: integer
% - Parameter jMode: non-negative integer
% - Returns index: a non-negative integer
arguments (Input)
    self WVTransform {mustBeNonempty}
    kMode (:,1) double {mustBeInteger}
    lMode (:,1) double {mustBeInteger}
    jMode (:,1) double {mustBeInteger,mustBeNonnegative}
end
arguments (Output)
    bool (:,1) logical {mustBeMember(bool,[0 1])}
end
bool = self.isValidPrimaryModeNumber(kMode,lMode,jMode) | self.isValidConjugateModeNumber(kMode,lMode,jMode);
end