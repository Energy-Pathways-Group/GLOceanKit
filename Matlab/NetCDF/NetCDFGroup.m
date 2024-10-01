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
    end

    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % TODO: complete initialization process, pass parentID?

        function self = NetCDFGroup(name,id)
            % NetCDFGroup
            %
            % - Topic: Initializing
            % - Declaration: ncfile = NetCDFFile(path,options)
            % - Parameter path: path to write file
            % - Parameter shouldOverwriteExisting: (optional) boolean indicating whether or not to overwrite an existing file at the path. Default 0.
            % - Returns: a new NetCDFFile instance
            arguments
                name char {mustBeNonempty}
                id double {mustBeNonempty}
            end
            self.name = name;
            self.id = id;

            self.dimensions = NetCDFDimension.empty(0,0);
            self.variables = NetCDFVariable.empty(0,0);
            self.attributes = containers.Map;
            self.complexVariables = NetCDFComplexVariable.empty(0,0);

            self.dimensionIDMap = configureDictionary('double','NetCDFDimension');
            self.dimensionNameMap = configureDictionary('string','NetCDFDimension');

            self.variableIDMap = configureDictionary('double','NetCDFVariable');
            self.variableNameMap = configureDictionary('string','NetCDFVariable');

            if isfield(options,'id')
                self.initializeGroupFromFile(group,options.id);
            else
                requiredFields = {'name'};
                for iField=1:length(requiredFields)
                    if ~isfield(options,requiredFields{iField})
                        error('You must specify %s to initial a new variable.',requiredFields{iField});
                    end
                end
                self.initGroup(group,options.name);
            end
        end

        function initializeGroupFromFile(self, parentGroup)
            self.parentGroup = parentGroup;

            % Fetch all dimensions and convert to NetCDFDimension objects
            dimensionIDs = netcdf.inqDimIDs(self.id);
            for iDim=1:length(dimensionIDs)
                self.addDimensionPrimitive(NetCDFDimension(self,id=dimensionIDs(iDim)));
            end
            self.dimensions = reshape(self.dimensions,[],1);

            % Fetch all variables and convert to NetCDFVariable objects
            variableIDs = netcdf.inqVarIDs(self.id);
            for iVar=1:length(variableIDs)
                self.addVariablePrimitive(NetCDFVariable(self,id=variableIDs(iVar)));
            end
            self.variables = reshape(self.variables,[],1);

            % Grab all the global attributes
            [~,~,ngatts,~] = netcdf.inq(self.id);
            for iAtt=0:(ngatts-1)
                gattname = netcdf.inqAttName(self.id,netcdf.getConstant('NC_GLOBAL'),iAtt);
                self.attributes(gattname) = netcdf.getAtt(self.id,netcdf.getConstant('NC_GLOBAL'),gattname);
            end

            % Now check to see if any variables form a complex variable
            self.complexVariables = NetCDFComplexVariable.complexVariablesFromVariables(self.variables);

            childGroupIDs = netcdf.inqGrps(self.id);

        end

        function initGroup(parentGroup, name)

        end

        function addDimension(self,name,options)
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
        end

        function addDimensionPrimitive(self,dimension)
            arguments
                self (1,1) NetCDFGroup {mustBeNonempty}
                dimension (1,1) NetCDFDimension {mustBeNonempty}
            end
            self.dimensions(end+1) = dimension;
            self.dimensionIDMap(dimension.varID) = dimension;
            self.dimensionNameMap(dimension.name) = dimension;
        end

        function addVariablePrimitive(self,variable)
            arguments
                self (1,1) NetCDFGroup {mustBeNonempty}
                variable (1,1) NetCDFVariable {mustBeNonempty}
            end
            self.variables(end+1) = variable;
            self.variableIDMap(variable.varID) = variable;
            self.variableNameMap(variable.name) = variable;
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
    end
end