function names = classRequiredPropertyNames()
% return a cell array of property names required by the class
%
% This function returns an array of property names required to be written
% by the class, in order to restore its state.
%
% - Topic: Developer
% - Declaration:  names = classRequiredPropertyNames()
% - Returns names: array strings
arguments (Output)
    names cell
end
names = {};
end