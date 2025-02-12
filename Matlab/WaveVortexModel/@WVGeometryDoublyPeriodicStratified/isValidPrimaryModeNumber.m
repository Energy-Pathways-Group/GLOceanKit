function bool = isValidPrimaryModeNumber(self,kMode,lMode,jMode)
% returns a boolean indicating whether (k,l,j) is a valid primary (non-conjugate) mode number
%
% returns a boolean indicating whether (k,l,j) is a valid
% non-conjugate mode number according to how the property
% conjugateDimension is set.
%
% - Topic: Index Gymnastics
% - Declaration: index = isValidPrimaryModeNumber(kMode,lMode,jMode)
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
klCheck = self.isValidPrimaryKLModeNumber(kMode,lMode);
jCheck = jMode >= 0 & jMode <= self.Nj;
bool = klCheck & jCheck;
end