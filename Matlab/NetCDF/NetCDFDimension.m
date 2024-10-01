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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = NetCDFDimension(group,options)
            arguments
                group (1,1) NetCDFGroup
                options.name string
                options.nPoints (1,1) double {mustBeNonnegative,mustBeInteger}
                options.id (1,1) double
            end
            if isfield(options,'id')
                self.initDimensionFromFile(group,options.id);
            else
                requiredFields = {'name','nPoints'};
                for iField=1:length(requiredFields)
                    if ~isfield(options,requiredFields{iField})
                        error('You must specify %s to initial a new dimension.',requiredFields{iField});
                    end
                end
                self.initDimension(group,options.name,options.nPoints);
            end
        end

        function initDimensionFromFile(self,group,id)
            % initialize an existing dimension from file
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

        function initDimension(self,group,name,nPoints)
            % initialize a new dimension
            %
            % To create a new dimension, you must specify the group, the
            % name, and the length, where the length is either a finite
            % value or inf, if the variable is to be mutable.
            arguments
                self (1,1) NetCDFDimension {mustBeNonempty}
                group (1,1) NetCDFGroup {mustBeNonempty} 
                name string {mustBeNonempty}
                nPoints (1,1) {mustBeNonnegative}
            end

            self.group = group;
            self.name = name;
   
            if isinf(nPoints)
                self.isMutable = 1;
                self.nPoints = 0;
                self.id = netcdf.defDim(self.group.id, name, netcdf.getConstant('NC_UNLIMITED'));
            else
                self.nPoints = nPoints;
                self.id = netcdf.defDim(self.group.id, name, self.nPoints);
            end 
        end

    end
end