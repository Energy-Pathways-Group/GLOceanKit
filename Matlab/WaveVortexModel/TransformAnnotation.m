classdef TransformAnnotation < handle
    properties
        name
        description
        detailedDescription
    end

    methods
        function self = TransformAnnotation(name,description,options)
            arguments
                name char {mustBeNonempty}
                description char {mustBeNonempty}
                options.detailedDescription char = ''
            end
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            self.name = name;
            self.description = description;
            self.loadDetailedDescriptionIfAvailable;
            if ~isempty(options.detailedDescription)
                if ~isempty(self.detailedDescription)
                    warning('Founded a detailedDescription md file for %s, but one was also set in code!',self.name);
                end
                self.detailedDescription = options.detailedDescription;
            end
        end

        function loadDetailedDescriptionIfAvailable(self)
            % dir(strcat(fileparts(mfilename('fullpath')),'/**/x.md'))
            % But this won't work from within the TransformDimension definition. Need
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