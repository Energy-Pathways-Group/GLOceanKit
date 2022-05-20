classdef WaveVortexVariable < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        name
        dimensions
        units
        description
        isComplex = 0 % does it have a non-zero imaginary part?
        isDimension = 0
        isVariableWithLinearTimeStep = 0
        isVariableWithNonlinearTimeStep = 0
    end

    methods
        function self = WaveVortexVariable(name,dimensions,units,description)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            if ~iscell(dimensions)
                if isempty(dimensions)
                    dimensions = {};
                else
                    dimensions = {dimensions};
                end
            end
            self.name = name;
            self.dimensions = dimensions;
            self.units = units;
            self.description = description;
        end

        function set.isVariableWithLinearTimeStep(self,value)
            self.isVariableWithLinearTimeStep = value;
            if value == 1
                self.isVariableWithNonlinearTimeStep = 1;
            end
        end
    end
end