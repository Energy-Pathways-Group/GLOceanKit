classdef WVPropertyAnnotation < WVAnnotation
    %Describes a property of the WVTransform
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
    %
    % - Declaration: classdef WVPropertyAnnotation < [WVAnnotation](/classes/wvannotation/)
    properties (GetAccess=public, SetAccess=private)
        % ordered cell array with the names of the dimensions
        % 
        % If the property has no dimensions, and empty cell array should be
        % passed. The dimension names must correspond to existing
        % dimensions.
        % - Topic: Properties
        dimensions

        % units of the dimension
        %
        % All units should be abbreviated SI units, e.g., 'm', or 'rad'.
        % - Topic: Properties
        units

        % boolean indicating whether or not the property may have an imaginary part
        %
        % This information is used when allocating space in a NetCDF file.
        % - Topic: Properties
        isComplex = 0
    end

    methods
        function self = WVPropertyAnnotation(name,dimensions,units,description,options)
            % create a new instance of WVPropertyAnnotation
            %
            % If a markdown file of the same name is in the same directory
            % or child directory, it will be loaded as the detailed
            % description upon initialization.
            %
            % - Topic: Initialization
            % - Declaration: propAnnotation = WVPropertyAnnotation(name,dimensions,units,description,options)
            % - Parameter name: name of the property
            % - Parameter dimensions: ordered list of the dimensions, or empty cell array
            % - Parameter units: abbreviated SI units of the property
            % - Parameter description: short description of the property
            % - Parameter isComplex: (optional) indicates whether the property has an imaginary part (default 0)
            % - Parameter detailedDescription: (optional) detailed description of the property
            % - Returns propAnnotation: a new instance of WVPropertyAnnotation
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