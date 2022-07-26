classdef WVVariableAnnotation < WVAnnotation
    % A variable describing the state of ocean in the WaveVortexTransform
    %
    % 
    %
    % Note that as a subclass of WVAnnotation, this class looks for
    % a file (name).md in the directory where it is defined another other
    % subdirectories. This file is then read-in to the detailed description
    % that is used on the website.

    properties
        dimensions
        units
        modelOp
        isComplex = 0 % does it have a non-zero imaginary part?
        isVariableWithLinearTimeStep = 1
        isVariableWithNonlinearTimeStep = 1
    end

    methods
        function self = WVVariableAnnotation(name,dimensions,units,description,options)
            arguments
                name char {mustBeNonempty}
                dimensions
                units char {mustBeNonempty}
                description char {mustBeNonempty}
                options.isComplex double {mustBeMember(options.isComplex,[0 1])} = 0
                options.detailedDescription char = ''
            end
            self@WVAnnotation(name,description,detailedDescription=options.detailedDescription);
            if ~iscell(dimensions)
                if isempty(dimensions)
                    dimensions = {};
                else
                    dimensions = {dimensions};
                end
            end
            self.dimensions = dimensions;
            self.units = units;
            self.isComplex = options.isComplex;
            if isempty(self.detailedDescription)
                self.detailedDescription = "- Topic: State Variables";
            end
        end

    end
end