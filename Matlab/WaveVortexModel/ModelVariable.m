classdef ModelVariable < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        name
        dimensions
        units
        description
        modelOp
        isComplex = 0 % does it have a non-zero imaginary part?
        isVariableWithLinearTimeStep = 1
        isVariableWithNonlinearTimeStep = 1
    end

    methods
        function self = ModelVariable(name,dimensions,units,description)
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

%         function set.isVariableWithLinearTimeStep(self,value)
%             self.isVariableWithLinearTimeStep = value;
%             if value == 1
%                 self.isVariableWithNonlinearTimeStep = 1;
%             end
%         end
    end
end