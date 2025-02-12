function names = spectralDimensionNames(self)
% return a cell array of the spectral dimension names
%
% This function returns an array of dimension names
%
% - Topic: Developer
% - Declaration:  names = spatialDimensionNames()
% - Returns names: array strings
arguments (Input)
    self WVTransform
end
arguments (Output)
    names cell
end
names = feval(strcat(class(self),'.classDefinedSpectralDimensionNames'));
end