classdef WVDimensionAnnotation < WVAnnotation
    %Describes a coordinate dimension of the WVTransform
    %
    % In addition to adding a name, description and detailed description of
    % a given dimension, you also specify its units. These annotations
    % are used for both online documentation and for writing to NetCDF
    % files.
    %
    % Note that as a subclass of WVAnnotation, this class looks for
    % a file (name).md in the directory where it is defined another other
    % subdirectories. This file is then read-in to the detailed description
    % that is used on the website.
    %
    % - Declaration: classdef WVDimensionAnnotation < [WVAnnotation](/classes/wvannotation/)
    properties
        % units of the dimension
        % 
        % All units should be abbreviated SI units, e.g., 'm', or 'rad'.
        % - Topic: Properties
        units
    end

    methods
        function self = WVDimensionAnnotation(name,units,description)
            % create a new instance of WVDimensionAnnotation
            %
            % If a markdown file of the same name is in the same directory
            % or child directory, it will be loaded as the detailed
            % description upon initialization.
            %
            % - Topic: Initialization
            % - Declaration: dimAnnotation = WVDimensionAnnotation(name,description,options)
            % - Parameter name: name of dimension
            % - Parameter units: abbreviated SI units of the dimension
            % - Parameter description: short description of the dimension
            % - Returns dimAnnotation: a new instance of WVDimensionAnnotation
            arguments
                name char {mustBeNonempty}
                units char {mustBeNonempty}
                description char {mustBeNonempty}
            end
            self@WVAnnotation(name,description);
            self.units = units;
        end
    end
end