classdef NetCDFDimension < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        group
        id
        name
        nPoints
        isMutable = 0
    end

    methods
        function self = NetCDFDimension(group,options)
            arguments
                group (1,1) NetCDFGroup
                options.name string
                options.nPoints (1,1) double {mustBeNonnegative,mustBeInteger}
                options.isMutable (1,1) double {mustBeMember(options.isMutable,[0 1])} = 0
                options.id (1,1) double
            end
            if isfield(options,'id')
                self.initDimensionFromFile(group,options.id);
            else
                requiredFields = {'name','nPoints','isMutable'};
                for iField=1:length(requiredFields)
                    if ~isfield(options,requiredFields{iField})
                        error('You must specify %s to initial a new dimension.',requiredFields{iField});
                    end
                end
                self.initDimension(group,options.name,options.dimensions,options.attributes,options.type);
            end
        end

        function initDimensionFromFile(self,group,id)
            arguments
                self (1,1) NetCDFDimension {mustBeNonempty}
                group (1,1) NetCDFGroup {mustBeNonempty}
                id (1,1) double {mustBeNonnegative,mustBeInteger}
            end
            self.group = group;
            self.id = id;

            unlimitedDimIDs = netcdf.inqUnlimDims(self.group.id);
            [self.name,self.nPoints] = netcdf.inqDim(self.group.id,self.id);
            if any(ismember(self.id,unlimitedDimIDs))
                self.isMutable = 1;
            end
        end

        function initDimension(self,group,name,length)
            arguments
                self (1,1) NetCDFDimension {mustBeNonempty}
                group (1,1) NetCDFGroup {mustBeNonempty} 
                name string {mustBeNonempty}
                length (1,1) {mustBeNonnegative}
            end

            self.group = group;
            self.name = name;
            self.nPoints = length;

            if isinf(length)
                length = netcdf.getConstant('NC_UNLIMITED');
                self.isMutable = 1;
            end

            self.id = netcdf.defDim(self.group.id, name, length);
        end
    end
end