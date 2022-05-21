classdef ModelDimension < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        name
        nPoints
        units
        description
    end

    methods
        function self = ModelDimension(name,nPoints,units,description)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            self.name = name;
            self.nPoints = nPoints;
            self.units = units;
            self.description = description;
        end
    end
end