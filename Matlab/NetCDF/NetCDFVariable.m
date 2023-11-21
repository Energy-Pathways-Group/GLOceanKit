classdef NetCDFVariable < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        varID
        name
        dimensions
        attributes
        type
        isComplex=0
    end

    methods
        function self = NetCDFVariable(name,dimensions,attributes,type,varID)
            arguments
                name string
                dimensions (:,1) NetCDFDimension
                attributes containers.Map
                type string
                varID (1,1) double
            end
            self.name = name;
            self.dimensions = dimensions;
            self.attributes = attributes;
            self.varID = varID;
            self.type = type;
        end

    end
end