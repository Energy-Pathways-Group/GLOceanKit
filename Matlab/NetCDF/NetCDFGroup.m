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
        realVariableIDMap
        realVariableNameMap
        complexVariableNameMap
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
            self.variables = NetCDFRealVariable.empty(0,0);
            self.groups = NetCDFGroup.empty(0,0);
            self.attributes = containers.Map;
            self.complexVariables = NetCDFComplexVariable.empty(0,0);

            self.dimensionIDMap = configureDictionary('double','NetCDFDimension');
            self.dimensionNameMap = configureDictionary('string','NetCDFDimension');

            self.realVariableIDMap = configureDictionary('double','NetCDFRealVariable');
            self.realVariableNameMap = configureDictionary('string','NetCDFRealVariable');
            self.complexVariableNameMap = configureDictionary('string','NetCDFComplexVariable');

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
            % initialize an existing group from file
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
                self.addDimensionPrimitive(self.parentGroup.dimensions);
            end
            dimensionIDs = netcdf.inqDimIDs(self.id);
            for iDim=1:length(dimensionIDs)
                self.addDimensionPrimitive(NetCDFDimension(self,id=dimensionIDs(iDim)));
            end

            % Fetch all variables and convert to NetCDFRealVariable objects
            variableIDs = netcdf.inqVarIDs(self.id);
            for iVar=1:length(variableIDs)
                self.addRealVariablePrimitive(NetCDFRealVariable(self,id=variableIDs(iVar)));
            end

            % Grab all the global attributes
            [~,~,ngatts,~] = netcdf.inq(self.id);
            for iAtt=0:(ngatts-1)
                gattname = netcdf.inqAttName(self.id,netcdf.getConstant('NC_GLOBAL'),iAtt);
                self.attributes(gattname) = netcdf.getAtt(self.id,netcdf.getConstant('NC_GLOBAL'),gattname);
            end

            % Now check to see if any variables form a complex variable
            self.addComplexVariablePrimitive(NetCDFComplexVariable.complexVariablesFromVariables(self.variables));
            if ~isempty(self.complexVariables)
                self.removeVariablePrimitive([self.complexVariables.realp self.complexVariables.imagp]);
            end

            groupIDs = netcdf.inqGrps(self.id);
            for iGrp=1:length(groupIDs)
                self.addGroupPrimitive(NetCDFGroup(parentGroup=self,id=groupIDs(iGrp)));
            end
        end

        function initGroup(self,parentGroup, name)
            % initialize a new group
            arguments
                self (1,1) NetCDFGroup {mustBeNonempty} 
                parentGroup (1,1) NetCDFGroup {mustBeNonempty} 
                name string {mustBeText}
            end
            self.parentGroup = parentGroup;
            self.name = name;
            self.id = netcdf.defGrp(self.parentGroup.id,self.name);
            self.addDimensionPrimitive(self.parentGroup.dimensions);
        end

        function grp = addGroup(self, name)
            arguments
                self (1,1) NetCDFGroup {mustBeNonempty} 
                name string {mustBeText}
            end
            if isKey(self.groupNameMap,name)
                error('A group with that name already exists.');
            end
            grp = NetCDFGroup(parentGroup=self,name=name);
            self.addGroupPrimitive(grp);
        end

        function [dim,var] = addDimension(self,name,value,options)
            % add a new dimension to the group
            %
            % This function adds both a dimension and an associated
            % coordinate variable.
            %
            % You may initialize directly with data, e.g.,
            %
            % ```matlab
            % ncfile.addDimension('x',0:9);
            % ```
            %
            % which will immediately add the dimension and variable data to
            % file, OR you may initialize the dimension and an empty
            % variable, e.g.,
            %
            % ```matlab
            % ncfile.addDimension('t',length=10,type='double');
            % ```
            %
            % which may be useful when initializing an unlimited
            % dimension. In this case you may also want to specify whether
            % or not the variable will be complex valued (it will default
            % to assuming the variable will be real valued).
            arguments
                self (1,1) NetCDFGroup {mustBeNonempty}
                name {mustBeText}
                value = []
                options.length = 0
                options.type char = []
                options.attributes containers.Map = containers.Map(KeyType='char',ValueType='any');
            end
            if isKey(self.dimensionNameMap,name) || isKey(self.realVariableNameMap,name)
                error('A dimension with that name already exists.');
            end
            if ~( (options.length > 0 && ~isempty(options.type)) || ~isempty(value))
                error('You must specify either the value of the data, or the length and type.');
            end

            if ~isempty(value)
                n = length(value);
            else
                n = options.length;
            end
            dim = NetCDFDimension(self,name=name,nPoints=n);
            self.addDimensionPrimitive(dim)
            % TODO: This needs to pass this down to child groups

            if ~isempty(value)
                options.type = class(value);
                options.isComplex = ~isreal(value);
            end
            var = self.addVariable(name,{dim.name},value,type=options.type,isComplex=0,attributes=options.attributes);           
        end

        function var = addVariable(self,name,dimNames,value,options)
            % add a new (real or complex) variable to the file
            %
            % Immediately initializes and writes the given variable data to
            % file, e.g.,
            %
            % ```matlab
            % ncfile.addVariable('fluid-tracer', {'x','y','t'}, myVariableData);
            % ```
            %
            % or intializes an variable without setting the data,
            % ```matlab
            % ncfile.addVariable('fluid-tracer', {'x','y','t'}, myVariableData);
            % ```
            %
            % - Topic: Working with variables
            % - Declaration: variable = addVariable(name,data,dimNames,properties,ncType)
            % - Parameter name: name of the variable (a string)
            % - Parameter data: variable data
            % - Parameter dimNames: cell array containing the dimension names
            % - Parameter properties: (optional) `containers.Map`
            % - Parameter properties: ncType
            % - Returns variable: a NetCDFVariable object
            arguments
                self (1,1) NetCDFGroup {mustBeNonempty}
                name {mustBeText}
                dimNames = {}
                value {mustBeNumericOrLogical} = [] 
                options.type char = []
                options.isComplex logical {mustBeMember(options.isComplex,[0 1])}
                options.attributes containers.Map = containers.Map(KeyType='char',ValueType='any');
            end
            if isKey(self.realVariableNameMap,name) || isKey(self.complexVariableNameMap,name)
                error('A variable with that name already exists.');
            end
            if ~( (~isempty(options.type) && ~isempty(options.isComplex)) || ~isempty(value))
                error('You must specify either the value of the data, or the ncType and isComplex.');
            end

            if ~isempty(value)
                options.type = class(value);
                options.isComplex = ~isreal(value);
            end

            if options.isComplex==1
                var = NetCDFComplexVariable(group=self,name=name,dimensions=self.dimensionWithName(dimNames),attributes=options.attributes,type=options.type);
                self.addComplexVariablePrimitive(var);
            else
                var = NetCDFRealVariable(self,name=name,dimensions=self.dimensionWithName(dimNames),attributes=options.attributes,type=options.type);
                self.addRealVariablePrimitive(var);
            end
            if ~isempty(value)
                var.value = value;
            end

        end

        function varargout = variable(self,variableNames)
            % read variables from file
            %
            % Pass a list of variables to read and the data will be
            % returned in the same order.
            %
            % ```matlab
            % [x,y] = ncfile.variable('x','y');
            % ```
            %
            % - Topic: Working with variables
            % - Declaration: varargout = variable(variableNames)
            % - Parameter variableNames: (repeating) list of variable names
            % - Returns varargout: (repeating) list of variable data
            arguments
                self NetCDFFile {mustBeNonempty}
            end
            arguments (Repeating)
                variableNames char
            end
            varargout = cell(size(variableNames));
            for iArg=1:length(variableNames)
                if isKey(self.complexVariableNameMap,variableNames{iArg})
                    varargout{iArg} = self.complexVariableNameMap(variableNames{iArg}).value;
                elseif isKey(self.realVariableNameMap,variableNames{iArg})
                    varargout{iArg} = self.realVariableNameMap(variableNames{iArg}).value;
                else
                    error('Unable to find a variable with the name %s',variableNames{iArg});
                end
            end
        end

        function varargout = variableAtIndexAlongDimension(self,dimName,index,variableNames)
            % read variables from file at a particular index (e.g., time)
            %
            % Pass a list of variables to read and the data will be
            % returned in the same order.
            %
            % ```matlab
            % [u,v] = ncfile.readVariablesAtIndexAlongDimension('t',100,'u','v');
            % ```
            %
            % - Topic: Working with variables
            % - Declaration: varargout = readVariables(variableNames)
            % - Parameter dimName: name of the dimension, character string
            % - Parameter index: index at which to read the data, positive integer
            % - Parameter variableNames: (repeating) list of variable names
            % - Returns varargout: (repeating) list of variable data
            arguments
                self NetCDFFile {mustBeNonempty}
                dimName char {mustBeNonempty}
                index  (1,1) double {mustBePositive} = 1
            end
            arguments (Repeating)
                variableNames char
            end
            varargout = cell(size(variableNames));
            for iArg=1:length(variableNames)
                if isKey(self.complexVariableNameMap,variableNames{iArg})
                    varargout{iArg} = self.complexVariableNameMap(variableNames{iArg}).valueAlongDimensionAtIndex(dimName,index);
                elseif isKey(self.realVariableNameMap,variableNames{iArg})
                    varargout{iArg} = self.realVariableNameMap(variableNames{iArg}).valueAlongDimensionAtIndex(dimName,index);
                else
                    error('Unable to find a variable with the name %s',variableNames{iArg});
                end
            end
        end

        function addDimensionPrimitive(self,dimension)
            arguments
                self (1,1) NetCDFGroup {mustBeNonempty}
                dimension (:,1) NetCDFDimension {mustBeNonempty}
            end
            for iDim=1:length(dimension)
                self.dimensions(end+1) = dimension(iDim);
                self.dimensionIDMap(dimension.id) = dimension(iDim);
                self.dimensionNameMap(dimension.name) = dimension(iDim);
            end
            self.dimensions = reshape(self.dimensions,[],1);
        end

        function addRealVariablePrimitive(self,variable)
            arguments
                self (1,1) NetCDFGroup {mustBeNonempty}
                variable (:,1) NetCDFRealVariable {mustBeNonempty}
            end
            for iVar=1:length(variable)
                self.variables(end+1) = variable(iVar);
                self.realVariableIDMap(variable.id) = variable(iVar);
                self.realVariableNameMap(variable.name) = variable(iVar);
            end
            self.variables = reshape(self.variables,[],1);
        end

        function removeVariablePrimitive(self,variable)
            arguments
                self (1,1) NetCDFGroup {mustBeNonempty}
                variable (:,1) NetCDFRealVariable
            end
            for iVar=1:length(variable)
                self.variables(self.variables==variable(iVar)) = [];
                self.realVariableIDMap(variable(iVar).id) = [];
                self.realVariableNameMap(variable(iVar).name) = [];
            end
        end

        function addComplexVariablePrimitive(self,variable)
            arguments
                self (1,1) NetCDFGroup {mustBeNonempty}
                variable (:,1) NetCDFComplexVariable
            end
            for iVar=1:length(variable)
                self.complexVariables(end+1) = variable(iVar);
                self.complexVariableNameMap(variable(iVar).name) = variable(iVar);
            end
            self.complexVariables = reshape(self.complexVariables,[],1);
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


        % key-value Map to retrieve a NetCDFRealVariable object by name
        % - Topic: Working with variables
        function v = variableWithName(self,name)
            v = self.realVariableNameMap(name);
        end


        % key-value Map to retrieve a NetCDFComplexVariable object by name
        % - Topic: Working with variables
        function v = complexVariableWithName(self,name)

        end

        % key-value Map to retrieve a NetCDFGroup object by name
        % - Topic: Working with groups
        function g = groupWithName(self,name)

        end

        function dump(self,options)
            arguments
                self NetCDFGroup
                options.indentLevel = 0
            end
            indent0 = sprintf(repmat('   ',1,options.indentLevel));
            indent1 = sprintf(repmat('   ',1,options.indentLevel+1));
            indent2 = sprintf(repmat('   ',1,options.indentLevel+2));
            indent3 = sprintf(repmat('   ',1,options.indentLevel+3));

            if isempty(self.parentGroup) && isa(self,'NetCDFFile')
                [~,filename,~] = fileparts(self.path);
                fprintf('%snetcdf: %s {\n',indent0,filename);
            else
                fprintf('\n%sgroup: %s {\n',indent0,self.name);
            end

            if ~isempty(self.dimensions)
                fprintf('%sdimensions: \n',indent1);
                for iDim=1:length(self.dimensions)
                    dimension = self.dimensions(iDim);
                    fprintf('%s%s = %d\n',indent2,dimension.name,dimension.nPoints);
                end
            end

            if ~isempty(self.variables)
                fprintf('\n%svariables: \n',indent1);
                for iVar=1:length(self.variables)
                    variable = self.variables(iVar);
                    fprintf('%s%s %s(',indent2,variable.type,variable.name);
                    for iDim=1:length(variable.dimensions)
                        if iDim==length(variable.dimensions)
                            fprintf('%s',variable.dimensions(iDim).name);
                        else
                            fprintf('%s,',variable.dimensions(iDim).name);
                        end
                    end
                    fprintf(')\n');
                    for attName = string(variable.attributes.keys)
                        if isa(variable.attributes(attName),'char') || isa(variable.attributes(attName),'string')
                            fprintf('%s%s = \"%s\"\n',indent3,attName,variable.attributes(attName));
                        elseif isa(variable.attributes(attName),'double') || isa(variable.attributes(attName),'single')
                            fprintf('%s%s = %f\n',indent3,attName,variable.attributes(attName));
                        else
                            fprintf('%s%s = %d\n',indent3,attName,variable.attributes(attName));
                        end
                    end
                end
            end

            if ~isempty(self.complexVariables)
                % fprintf('\n%scomplex variables: \n',indent1);
                for iVar=1:length(self.complexVariables)
                    variable = self.complexVariables(iVar);
                    fprintf('%scomplex %s %s(',indent2,variable.type,variable.name);
                    for iDim=1:length(variable.dimensions)
                        if iDim==length(variable.dimensions)
                            fprintf('%s',variable.dimensions(iDim).name);
                        else
                            fprintf('%s,',variable.dimensions(iDim).name);
                        end
                    end
                    fprintf(')\n');
                    for attName = string(variable.attributes.keys)
                        if isa(variable.attributes(attName),'char') || isa(variable.attributes(attName),'string')
                            fprintf('%s%s = \"%s\"\n',indent3,attName,variable.attributes(attName));
                        elseif isa(variable.attributes(attName),'double') || isa(variable.attributes(attName),'single')
                            fprintf('%s%s = %f\n',indent3,attName,variable.attributes(attName));
                        else
                            fprintf('%s%s = %d\n',indent3,attName,variable.attributes(attName));
                        end
                    end
                end
            end

            if ~isempty(self.groups)
                for iGroup=1:length(self.groups)
                    group = self.groups(iGroup);
                    group.dump(indentLevel=options.indentLevel+1);
                end
            end

            if self.attributes.Count > 0
                fprintf('\n%sglobal attributes: \n',indent1);
                for attName = string(self.attributes.keys)
                    if isa(self.attributes(attName),'char') || isa(self.attributes(attName),'string')
                        fprintf('%s%s = \"%s\"\n',indent2,attName,self.attributes(attName));
                    elseif isa(self.attributes(attName),'double') || isa(self.attributes(attName),'single')
                        fprintf('%s%s = %f\n',indent2,attName,self.attributes(attName));
                    else
                        fprintf('%s%s = %d\n',indent2,attName,self.attributes(attName));
                    end
                end
            end

            if isempty(self.parentGroup) && isa(self,'NetCDFFile')
                fprintf('%s}\n',indent0);
            else
                fprintf('%s} // group %s\n',indent0,self.name);
            end
        end
    end
end