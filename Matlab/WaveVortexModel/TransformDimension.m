classdef TransformDimension < TransformAnnotation
    %Describes a coordinate dimension of the WaveVortexTransform
    %
    % In addition to adding a name, description and detailed description of
    % a given dimension, you also specify tits  units. These annotations
    % are used for both online documentation and for writing to NetCDF
    % files.
    %
    % Note that as a subclass of TransformAnnotation, this class looks for
    % a file (name).md in the directory where it is defined another other
    % subdirectories. This file is then read-in to the detailed description
    % that is used on the website.
    properties
        units
    end

    methods
        function self = TransformDimension(name,units,description)
            arguments
                name char {mustBeNonempty}
                units char {mustBeNonempty}
                description char {mustBeNonempty}
            end
            self@TransformAnnotation(name,description);
            self.units = units;
        end
    end
end