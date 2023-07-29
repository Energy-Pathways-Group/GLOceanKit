classdef WVAnnotation < handle
    % annotates methods, properties, operations, and variables
    %
    % The purpose of this class is twofold. First, it lets us add
    % descriptions and detailedDescriptions for any method or property both
    % inline in code, and externally in a markdown file. Second, its
    % subclasses support annotating units and dimensions, which is not only
    % useful for help documentation, but is also necessary for adding those
    % variables to NetCDF files.
    %
    % This class looks for a detailedDescription in a .md file with the
    % same name. It will look up to one folder deep from the call site.
    %
    % - Declaration: classdef WVAnnotation < handle
    properties (GetAccess=public, SetAccess=private)
        % name of the method, property, or variable
        % 
        % The name should be an exact match to the method, property, or
        % variable that it is describing.
        % - Topic: Properties
        name

        % short description of the method, property, or variable
        % 
        % The short description is used in the table-of-contents of
        % documentation and as the long_name property in NetCDF files.
        % - Topic: Properties
        description
    end

    properties (GetAccess=public, SetAccess=public)
        % a detailed description of the method, property, or variable
        % 
        % The detailed description may be written in markdown and may also
        % use MathJax latex notation. The detailed description will be used
        % for the help documentation.
        %
        % The detailed description can be populated using a markdown
        % sidecar file. Specifically, when the WVAnnotation (or subclass)
        % is initialized, it will look in the directory and one directory
        % deep for a markdown file with the name of the annotation and the
        % file extension .md.
        % - Topic: Properties
        detailedDescription

        attributes
    end

    methods
        function self = WVAnnotation(name,description,options)
            % create a new instance of WVAnnotation
            %
            % Creates a new instance of WVAnnotation with a name,
            % description and optional detailed description.
            %
            % If a markdown file of the same name is in the same directory
            % or child directory, it will be loaded as the detailed
            % description upon initialization.
            %
            % - Topic: Initialization
            % - Declaration: wvAnnotation = WVAnnotation(name,description,options)
            % - Parameter name: name of the method, property, or variable
            % - Parameter description: short description of the method, property, or variable
            % - Parameter detailedDescription: (optional) a detailed description of the method, property, or variable
            % - Returns wvAnnotation: a new instance of WVAnnotation
            arguments
                name char {mustBeNonempty}
                description char {mustBeNonempty}
                options.detailedDescription char = ''
                options.attributes = containers.Map();
            end
            self.name = name;
            self.description = description;
            self.loadDetailedDescriptionIfAvailable;
            self.attributes = options.attributes;
            if ~isempty(options.detailedDescription)
%                 if ~isempty(self.detailedDescription)
%                     warning('Founded a detailedDescription md file for %s, but one was also set in code!',self.name);
%                 end
                self.detailedDescription = options.detailedDescription;
            end
        end
    end

    methods (Access=protected)
        function loadDetailedDescriptionIfAvailable(self)
            % dir(strcat(fileparts(mfilename('fullpath')),'/**/x.md'))
            % But this won't work from within the WVDimensionAnnotation definition. Need
            % to use dbstack to get the caller name.
            st = dbstack;
            if length(st) > 2
                files = dir(strcat(fileparts(which(st(3).file)),'/**/',self.name,'.md'));
                if ~isempty(files)
                    self.detailedDescription = fileread(strcat(files(1).folder,'/',files(1).name));
                end
            end
        end
    end
end