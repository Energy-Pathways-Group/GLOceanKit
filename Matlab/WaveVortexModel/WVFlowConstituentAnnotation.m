classdef WVFlowConstituentAnnotation < WVAnnotation
    %Describes a set of energetically orthogonal flow states
    %
    % 
    %
    % - Declaration: classdef WVFlowConstituentAnnotation < [WVAnnotation](/classes/wvannotation/)
    properties
        % name of the major constituent group
        % 
        % For example, internal gravity waves belong to the "waves" group.
        % - Topic: Properties
        groupName

        % abbreviated name
        % 
        % abreviated name, e.g., "igw" for internal gravity waves.
        % - Topic: Properties
        abbreviatedName

        
    end

    methods
        function self = WVFlowConstituentAnnotation(name,groupName,description)
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
                groupName char {mustBeNonempty}
                description char {mustBeNonempty}
            end
            self@WVAnnotation(name,description);
            self.groupName = groupName;
        end
    end
end