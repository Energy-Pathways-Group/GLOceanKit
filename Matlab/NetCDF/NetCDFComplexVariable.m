classdef NetCDFComplexVariable < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        name
        realVar
        imagVar
    end

    methods
        function self = NetCDFComplexVariable(name,realVar,imagVar)
            self.name = name;
            self.realVar = realVar;
            self.imagVar = imagVar;
        end

    end
end