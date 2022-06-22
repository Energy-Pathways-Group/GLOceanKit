classdef StateVariable < handle
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
        function self = StateVariable(name,dimensions,units,description,options)
            arguments
                name char {mustBeNonempty}
                dimensions
                units char {mustBeNonempty}
                description char {mustBeNonempty}
                options.isComplex double {mustBeMember(options.isComplex,[0 1])} = 0
            end
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
            self.isComplex = options.isComplex;
        end

%         function set.isVariableWithLinearTimeStep(self,value)
%             self.isVariableWithLinearTimeStep = value;
%             if value == 1
%                 self.isVariableWithNonlinearTimeStep = 1;
%             end
%         end
    end
end