classdef NetCDFGroup < handle
    properties
        % id of the group
        % - Topic: Accessing group properties
        id

        % name of the group
        % - Topic: Accessing group properties
        name

        parentGroup

        % array of NetCDFGroup objects
        % - Topic: Working with groups
        groups

        % array of NetCDFDimension objects
        %
        % An array of NetCDFDimension objects for each coordinate dimension
        % defined in the NetCDF file.
        %
        % Usage
        % ```matlab
        % dim = ncfile.dimensions(dimID+1); % get the dimension with dimID
        % ```
        % - Topic: Working with dimensions
        dimensions

        % array of NetCDFVariable objects
        % - Topic: Working with variables
        variables

        % key-value Map of global attributes
        %
        % A `containers.Map` type that contains the key-value pairs of all
        % global attributes in the NetCDF file. This is intended to be
        % *read only*. If you need to add a new attribute to file, use
        % [`addAttribute`](#addattribute).
        %
        % Usage
        % ```matlab
        % model = ncfile.attributes('model');
        % ```
        % - Topic: Working with global attributes
        attributes

        % array of NetCDFComplexVariable objects
        % - Topic: Working with variables
        complexVariables            
    end

    properties (Access=private)
        dimensionIDMap
        dimensionNameMap
        variableIDMap
        variableNameMap
        groupIDMap
        groupNameMap
    end

    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % TODO: complete initialization process, pass parentID?

        function self = NetCDFGroup(options)
            % NetCDFGroup
            %
            % - Topic: Initializing
            % - Declaration: ncfile = NetCDFFile(path,options)
            % - Parameter path: path to write file
            % - Parameter shouldOverwriteExisting: (optional) boolean indicating whether or not to overwrite an existing file at the path. Default 0.
            % - Returns: a new NetCDFFile instance
            arguments
                options.name char {mustBeNonempty}
                options.id double {mustBeNonempty}
                options.parentGroup = []
            end
            
            self.dimensions = NetCDFDimension.empty(0,0);
            self.variables = NetCDFVariable.empty(0,0);
            self.groups = NetCDFGroup.empty(0,0);
            self.attributes = containers.Map;
            self.complexVariables = NetCDFComplexVariable.empty(0,0);

            self.dimensionIDMap = configureDictionary('double','NetCDFDimension');
            self.dimensionNameMap = configureDictionary('string','NetCDFDimension');

            self.variableIDMap = configureDictionary('double','NetCDFVariable');
            self.variableNameMap = configureDictionary('string','NetCDFVariable');

            self.groupIDMap = configureDictionary('double','NetCDFGroup');
            self.groupNameMap = configureDictionary('string','NetCDFGroup');

            if isfield(options,'id')
                self.initializeGroupFromFile(options.parentGroup,options.id);
            else
                requiredFields = {'name'};
                for iField=1:length(requiredFields)
                    if ~isfield(options,requiredFields{iField})
                        error('You must specify %s to initialize a new group.',requiredFields{iField});
                    end
                end
                self.initGroup(options.parentGroup,options.name);
            end
        end

        function initializeGroupFromFile(self,parentGroup,id)
            arguments
                self (1,1) NetCDFGroup {mustBeNonempty}
                parentGroup
                id (1,1) double {mustBeNonnegative,mustBeInteger}
            end
            self.id = id;
            self.parentGroup = parentGroup;
            self.name = netcdf.inqGrpName(self.id);

            % Fetch all dimensions and convert to NetCDFDimension objects
            if ~isempty(self.parentGroup)
                self.dimensions = self.parentGroup.dimensions;
            end
            dimensionIDs = netcdf.inqDimIDs(self.id);
            for iDim=1:length(dimensionIDs)
                self.addDimensionPrimitive(NetCDFDimension(self,id=dimensionIDs(iDim)));
            end

            % Fetch all variables and convert to NetCDFVariable objects
            variableIDs = netcdf.inqVarIDs(self.id);
            for iVar=1:length(variableIDs)
                self.addVariablePrimitive(NetCDFVariable(self,id=variableIDs(iVar)));
            end

            % Grab all the global attributes
            [~,~,ngatts,~] = netcdf.inq(self.id);
            for iAtt=0:(ngatts-1)
                gattname = netcdf.inqAttName(self.id,netcdf.getConstant('NC_GLOBAL'),iAtt);
                self.attributes(gattname) = netcdf.getAtt(self.id,netcdf.getConstant('NC_GLOBAL'),gattname);
            end

            % Now check to see if any variables form a complex variable
            self.complexVariables = NetCDFComplexVariable.complexVariablesFromVariables(self.variables);

            groupIDs = netcdf.inqGrps(self.id);
            for iGrp=1:length(groupIDs)
                self.addGroupPrimitive(NetCDFGroup(parentGroup=self,id=groupIDs(iGrp)));
            end
        end

        function initGroup(parentGroup, name)

        end

        function addDimension(self,name,options)
            % This is how the user should add a new dimension to the group
            arguments
                self (1,1) NetCDFDimension {mustBeNonempty}
                name string {mustBeNonempty}
                options.length = 0
                options.value = []
                options.attributes (1,1) containers.Map = containers.Map()
                options.ncType char = []
            end
            if ~isempty(group.dimensionWithName(name)) || ~isempty(group.variableWithName(name))
                error('A dimension with that name already exists.');
            end
            if (options.length == 0 && isempty(options.value))
                error('You must specify either the value or the length of the data');
            end
            if (isempty(options.value) && isempty(options.ncType))
                error('If you do not specify the value during initialization, you must specify the ncType');
            end

            % First create the dimension
            if isinf(options.length)
                n = inf;
            else
                n = length(options.data);
            end
            dim = NetCDFDimension(self,name=name,nPoints=n);
            self.addDimensionPrimitive(dim)

            % Now create the associated coordinate variable
            if isempty(options.value)
                ncType = options.ncType;
            else
                ncType = NetCDFFile.netCDFTypeForData(options.value);
            end
            var = NetCDFVariable(self,name=name,dimensions={dim},attributes=options.attributes,type=ncType);
            self.addVariablePrimitive(var);
            if ~isempty(options.value)
                var.value = options.value;
            end

            % TODO: This needs to pass this down to child groups
        end

        function addDimensionPrimitive(self,dimension)
            arguments
                self (1,1) NetCDFGroup {mustBeNonempty}
                dimension (1,1) NetCDFDimension {mustBeNonempty}
            end
            self.dimensions(end+1) = dimension;
            self.dimensions = reshape(self.dimensions,[],1);
            self.dimensionIDMap(dimension.id) = dimension;
            self.dimensionNameMap(dimension.name) = dimension;
        end

        function addVariablePrimitive(self,variable)
            arguments
                self (1,1) NetCDFGroup {mustBeNonempty}
                variable (1,1) NetCDFVariable {mustBeNonempty}
            end
            self.variables(end+1) = variable;
            self.variables = reshape(self.variables,[],1);
            self.variableIDMap(variable.id) = variable;
            self.variableNameMap(variable.name) = variable;
        end

        function addGroupPrimitive(self,group)
            arguments
                self (1,1) NetCDFGroup {mustBeNonempty}
                group (1,1) NetCDFGroup {mustBeNonempty}
            end
            self.groups(end+1) = group;
            self.groups = reshape(self.groups,[],1);
            self.groupIDMap(group.id) = group;
            self.groupNameMap(group.name) = group;
        end


        function variable = initVariable(self,name,dimNames,attributes,ncType)
            if isKey(self.variableNameMap,name)
                error('A variable with that name already exists.');
            end

            dims = self.dimensionWithName(dimNames);

        end

        % key-value Map to retrieve a NetCDFDimension object by name
        %
        % Usage
        % ```matlab
        % xDim = ncfile.dimensionWithName('x');
        % ```
        %
        % - Topic: Working with dimensions
        function dims = dimensionWithName(self,name)
            dims = self.dimensionNameMap(name);
        end

        function dims = dimensionWithID(self,dimids)
            % return the dimension IDs given the dimension names
            %
            % - Topic: Working with dimensions
            dims = self.dimensionIDMap(dimids);
        end


        % key-value Map to retrieve a NetCDFVariable object by name
        % - Topic: Working with variables
        function v = variableWithName(self,name)

        end


        % key-value Map to retrieve a NetCDFComplexVariable object by name
        % - Topic: Working with variables
        function v = complexVariableWithName(self,name)

        end

        % key-value Map to retrieve a NetCDFGroup object by name
        % - Topic: Working with groups
        function g = groupWithName(self,name)

        end

        function dump(self)
            for iVar=1:length(self.variables)
                variable = self.variables(iVar);
                fprintf('%s\t{',variable.name);
                for iDim=1:length(variable.dimensions)
                    if iDim==length(variable.dimensions)
                        fprintf('%s',variable.dimensions(iDim).name);
                    else
                        fprintf('%s,',variable.dimensions(iDim).name);
                    end
                end
                fprintf('}\n');
            end
        end
    end
end