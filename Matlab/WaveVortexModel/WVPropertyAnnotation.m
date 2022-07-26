classdef WVPropertyAnnotation < WVAnnotation
    %Describes a property of the WaveVortexTransform
    %
    % In addition to adding a name, description and detailed description of
    % a given property, you can also specify the properties dimensions,
    % its units, and whether it is a complex number or not. These
    % annotations are used for both online documentation and for writing to
    % NetCDF files.
    %
    % Note that as a subclass of WVAnnotation, this class looks for
    % a file (name).md in the directory where it is defined another other
    % subdirectories. This file is then read-in to the detailed description
    % that is used on the website.
    properties
        dimensions
        units
        isComplex = 0 % does it have a non-zero imaginary part?
    end

    methods
        function self = WVPropertyAnnotation(name,dimensions,units,description,options)
            arguments
                name char {mustBeNonempty}
                dimensions
                units char
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
        end
    end
end