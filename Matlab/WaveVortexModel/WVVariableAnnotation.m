classdef WVVariableAnnotation < WVAnnotation
    % Describes a variable computed from the WVTransform
    % 
    % In addition to adding a name, description and detailed description of
    % a given variable, you also specify its dimensions, units, and whether
    % or note it has an imaginary part. These annotations are used for both
    % online documentation and for writing to NetCDF files.
    %
    % Setting the two properties `isVariableWithLinearTimeStep` and
    % `isVariableWithNonlinearTimeStep` are important for determining
    % how the variable is cached, and when it is saved to a NetCDF file.
    %
    % Note that as a subclass of WVAnnotation, this class looks for
    % a file (name).md in the directory where it is defined another other
    % subdirectories. This file is then read-in to the detailed description
    % that is used on the website.
    %
    % - Declaration: classdef WVVariableAnnotation < [WVAnnotation](/classes/wvannotation)
    properties
        % ordered cell array with the names of the dimensions
        %
        % If the variable has no dimensions, and empty cell array should be
        % passed. The dimension names must correspond to existing
        % dimensions.
        % - Topic: Properties
        dimensions

        % units of the variable
        %
        % All units should be abbreviated SI units, e.g., 'm', or 'rad'.
        % - Topic: Properties
        units

        % WVOperation responsible for computing this variable
        %
        % This property will be automatically populated when the variable
        % annotation is passed to the WVOperation.
        % - Topic: Properties
        modelOp

        % boolean indicating whether or not the property may have an imaginary part
        %
        % This information is used when allocating space in a NetCDF file.
        % - Topic: Properties
        isComplex = 0

        % boolean indicating whether the variable changes value with a linear time step
        %
        % This information is used when caching variables and when writing
        % to NetCDF file.
        % - Topic: Properties
        isVariableWithLinearTimeStep = 1

        % boolean indicating whether the variable changes value with a non-linear time step
        %
        % This information is used when caching variables and when writing
        % to NetCDF file.
        % - Topic: Properties
        isVariableWithNonlinearTimeStep = 1
    end

    methods
        function self = WVVariableAnnotation(name,dimensions,units,description,options)
            % create a new instance of WVVariableAnnotation
            %
            % If a markdown file of the same name is in the same directory
            % or child directory, it will be loaded as the detailed
            % description upon initialization.
            %
            % - Topic: Initialization
            % - Declaration: variableAnnotation = WVVariableAnnotation(name,dimensions,units,description,options)
            % - Parameter name: name of the variable
            % - Parameter dimensions: ordered list of the dimensions, or empty cell array
            % - Parameter units: abbreviated SI units of the variable
            % - Parameter description: short description of the variable
            % - Parameter isComplex: (optional) indicates whether the variable has an imaginary part (default 0)
            % - Parameter detailedDescription: (optional) detailed description of the variable
            % - Returns variableAnnotation: a new instance of WVVariableAnnotation
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