classdef NetCDFDimension < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        dimID
        name
        nPoints
        isMutable = 0
    end

    methods
        function self = NetCDFDimension(name,nPoints,dimID)
            self.name = name;
            self.nPoints = nPoints;
            self.dimID = dimID;
        end

    end
end