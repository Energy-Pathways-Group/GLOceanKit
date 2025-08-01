classdef NetCDFDimension < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (WeakHandle)
        group NetCDFGroup
    end

    properties
        id
        name
        nPoints
        isMutable = 0
        namePath
    end

    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = NetCDFDimension(group,options)
            % create a NetCDFDimension
            %
            % All dimensions must have a group to which they belong. If you
            % initialize by passing an id, then the dimension will be
            % intialized from file. If you initialize by passing the name
            % and the number of points, then a new dimension will be
            % created and added to file.
            %
            % - Topic: Working with dimensions
            % - Declaration: self = NetCDFDimension(group,options)
            % - Parameter group: a NetCDFGroup instance
            % - Parameter name (optional): name of the dimension (string)
            % - Parameter nPoints (optional): 0 or more, inf if unlimited
            % - Parameter id (optional): id of an existing dimension
            % - Returns self: a NetCDFDimension object
            arguments
                group (1,1) NetCDFGroup
                options.name string
                options.nPoints (1,1) double {mustBeNonnegative}
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

        function updateLength(self)
            [~,self.nPoints] = netcdf.inqDim(self.group.id,self.id);
        end
        
        function val = get.namePath(self)
            if isequal(self.group.groupPath,"")
                val = self.name;
            else
                val = self.group.groupPath + "/" + self.name;
            end
        end

    end
end