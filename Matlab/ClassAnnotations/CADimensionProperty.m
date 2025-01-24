classdef CADimensionProperty < CAPropertyAnnotation
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
        % units of the dimension
        %
        % All units should be abbreviated SI units, e.g., 'm', or 'rad'.
        % - Topic: Properties
        units
    end

    methods
        function self = CADimensionProperty(name,units,description,options)
            % create a new instance of WVPropertyAnnotation
            %
            % If a markdown file of the same name is in the same directory
            % or child directory, it will be loaded as the detailed
            % description upon initialization.
            %
            % - Topic: Initialization
            % - Declaration: propAnnotation = WVPropertyAnnotation(name,dimensions,units,description,options)
            % - Parameter name: name of the property
            % - Parameter units: abbreviated SI units of the property
            % - Parameter description: short description of the property
            % - Parameter detailedDescription: (optional) detailed description of the property
            % - Returns propAnnotation: a new instance of WVPropertyAnnotation
            arguments
                name char {mustBeNonempty}
                units char
                description char {mustBeNonempty}
                options.detailedDescription char = ''
            end

            self@CAPropertyAnnotation(name,description,detailedDescription=options.detailedDescription);
            self.units = units;
        end
    end
end