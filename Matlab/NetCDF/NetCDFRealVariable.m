classdef NetCDFRealVariable < NetCDFVariable
    % Encapsulates variable data in a NetCDF file
    % 

    properties
        id
    end

    properties (Dependent, Hidden)
        value
    end

    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = NetCDFRealVariable(group,options)
            % initialize a NetCDFVariable either from file or scratch
            %
            % This function is not intended to be called directly by a
            % user. To create a new variable, the user should call
            % `addVariable' or `initVariable' on the group or file.
            %
            % To initialize a variable from file you must pass the group
            % and variable id, e.g.
            %
            % ```matlab
            % aVariable = NetCDFVariable(group,id=3);
            % ```
            %
            % To create a new variable, you must pass the group and the
            % other required parameters
            %
            % ```matlab
            % aVariable = NetCDFVariable(group,name=aName,attributes=someAttributes,type='double');
            % ```
            %
            % class(data)
            arguments
                group (1,1) NetCDFGroup
                options.name {mustBeText}
                options.dimensions (:,1) NetCDFDimension
                options.attributes containers.Map = containers.Map
                options.type string
                options.id (1,1) double
            end
            if isfield(options,'id')
                self.initVariableFromFile(group,options.id);
            else
                requiredFields = {'name','dimensions','type'};
                for iField=1:length(requiredFields)
                    if ~isfield(options,requiredFields{iField})
                        error('You must specify %s to initial a new variable.',requiredFields{iField});
                    end
                end
                self.initVariable(group,options.name,options.dimensions,options.type,options.attributes);
            end
        end

        function initVariableFromFile(self, group, id)
            % initialize a variable from file
            %
            % 
            arguments
                self (1,1) NetCDFRealVariable
                group (1,1) NetCDFGroup
                id (1,1) double
            end
            self.group = group;
            self.id = id;
            self.attributes = containers.Map;

            [self.name, xtype, dimids, numatts] = netcdf.inqVar(group.id,id);
            self.type = NetCDFVariable.matlabTypeForNetCDFType(NetCDFVariable.typeStringForTypeID(xtype));
            
            for iAtt=0:(numatts-1)
                gattname = netcdf.inqAttName(group.id,id,iAtt);
                self.attributes(gattname) = netcdf.getAtt(group.id,id,gattname);
            end
            self.dimensions = self.group.dimensionWithID(dimids);
        end

        function initVariable(self,group,name,dimensions,type,attributes)
            % initialize a real-valued variable
            %
            % The basic work flow is that you need to first,
            % - initVariable followed by,
            % - setVariable
            % Or you can just call
            % - addVariable
            % which will both initialize and then set. The reason to
            %
            % ```matlab
            % ncfile.initVariable('fluid-tracer', {'x','y','t'},'double',containers.Map({'isTracer'},{'1'}));
            % ncfile.setVariable('fluid-tracer',myVariable);
            % ```
            %
            % - Topic: Working with variables
            % - Declaration: variable = initVariable(name,dimNames,properties,ncType)
            % - Parameter name: name of the variable (a string)
            % - Parameter dimNames: cell array containing the dimension names
            % - Parameter properties: (optional) `containers.Map`
            % - Parameter properties: ncType
            % - Returns variable: a NetCDFVariable object
            arguments
                self (1,1) NetCDFRealVariable {mustBeNonempty}
                group (1,1) NetCDFGroup {mustBeNonempty}
                name {mustBeText}
                dimensions (:,1) NetCDFDimension
                type string {mustBeMember(type,["double" "single" "int64" "uint64" "int32" "uint32" "int16" "uint16" "int8" "uint8" "char" "string" "logical"])}
                attributes containers.Map = containers.Map(KeyType='char',ValueType='any')         
            end
            self.group = group;
            self.dimensions = dimensions;
            self.id = netcdf.defVar(group.id, name, NetCDFVariable.ncTypeForMatlabType(type), [self.dimensions.id]);
            self.name = name;
            self.type = type;

            isFunctionHandle = attributes.isKey(NetCDFVariable.GLNetCDFSchemaIsFunctionHandleTypeKey) && ~attributes(NetCDFVariable.GLNetCDFSchemaIsFunctionHandleTypeKey);
            if strcmp(type,"logical") %&& isFunctionHandle
                attributes(NetCDFVariable.GLNetCDFSchemaIsLogicalTypeKey) = uint8(1);
            end

            keyNames = keys(attributes);
            for iKey = 1:length(keyNames)
                netcdf.putAtt(group.id,self.id, keyNames{iKey}, attributes(keyNames{iKey}));
            end
            self.attributes = attributes;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Add attributes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function addAttribute(self,key,value)
            self.attributes(key) = value;
            netcdf.putAtt(self.group.id,self.id, key, value);
            % netcdf.sync(self.group.id);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Read/write data
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function value = get.value(self)
            arguments
                self NetCDFRealVariable
            end
            value = netcdf.getVar(self.group.id,self.id);
            if isvector(value) && isKey(self.attributes, self.GLNetCDFSchemaIsRowVectorKey)
                if self.attributes(self.GLNetCDFSchemaIsRowVectorKey) == 1
                    value=reshape(value,1,[]);
                else
                    value=reshape(value,[],1);
                end
            end
            if isKey(self.attributes,self.GLNetCDFSchemaIsLogicalTypeKey) && self.attributes(self.GLNetCDFSchemaIsLogicalTypeKey) == 1
                value = logical(value);
            end
            if isKey(self.attributes,self.GLNetCDFSchemaIsFunctionHandleTypeKey) && self.attributes(self.GLNetCDFSchemaIsFunctionHandleTypeKey) == 1
                binaryData = uint8(value);
                tmpfile = strcat(tempname,'.mat');
                fileID = fopen(tmpfile, 'w');
                fwrite(fileID,binaryData, '*uint8');
                fclose(fileID);
                matFile = load(tmpfile);
                value = matFile.(self.name);
                delete(tmpfile);
            end
        end

        function set.value(self,data)
            arguments
                self NetCDFRealVariable {mustBeNonempty}
                data {mustBeNonempty}
            end
            for iDim=1:length(self.dimensions)
                if (length(self.dimensions)==1 && length(data) ~= self.dimensions(iDim).nPoints) || (length(self.dimensions)~=1 && size(data,iDim) ~= self.dimensions(iDim).nPoints)
                    error('Incorrect dimension size: dimension %d of the data of %s is length %d, but the dimension %s has length %d.',iDim,self.name,size(data,iDim),self.dimensions(iDim).name,self.dimensions(iDim).nPoints);
                end
            end
            if isKey(self.attributes,self.GLNetCDFSchemaIsLogicalTypeKey) && self.attributes(self.GLNetCDFSchemaIsLogicalTypeKey) == 1
                data = uint8(data);
            end
            netcdf.putVar(self.group.id, self.id, data);
            if isvector(data)
                self.addAttribute(self.GLNetCDFSchemaIsRowVectorKey, uint8(isrow(data)));
            end
            % netcdf.sync(self.group.id);
        end

        function value = valueAlongDimensionAtIndex(self,dimensionName,index)
            arguments
                self NetCDFRealVariable {mustBeNonempty} 
                dimensionName char {mustBeNonempty}
                index {mustBeNonnegative}
            end
            dimension = self.group.dimensionWithName(dimensionName);
            start = zeros(1,length(self.dimensions));
            count = zeros(1,length(self.dimensions));
            for iDim=1:length(self.dimensions)
                if self.dimensions(iDim) == dimension
                    start(iDim) = index-1;
                    count(iDim) = 1;
                else
                    start(iDim) = 0;
                    count(iDim) = self.dimensions(iDim).nPoints;
                end
            end
            value = netcdf.getVar(self.group.id, self.id, start, count);
        end

        function setValueAlongDimensionAtIndex(self,data,dimensionName,index)
            % append new data to an existing variable
            %
            % concatenates data along a variable dimension (such as a time
            % dimension).
            %
            % ```matlab
            % variable.setValueAlongDimensionAtIndex(data,dimensionName,index);
            % ```
            %
            % - Topic: Working with variables
            % - Declaration: setValueAlongDimensionAtIndex(data,dimension,index)
            % - Parameter data: variable data
            % - Parameter dimension: the variable dimension along which to concatenate
            % - Parameter index: index at which to write data
            arguments
                self NetCDFRealVariable {mustBeNonempty} 
                data 
                dimensionName char {mustBeNonempty}
                index {mustBeNonnegative}
            end
            dimension = self.group.dimensionWithName(dimensionName);
            if any(dimension.nPoints < index) && dimension.isMutable ~= 1
                error('Cannot increase the length of a dimension that is not mutable');
            end

            start = zeros(1,length(self.dimensions));
            count = zeros(1,length(self.dimensions));
            for iDim=1:length(self.dimensions)
                if self.dimensions(iDim) == dimension
                    start(iDim) = index(1)-1;
                    count(iDim) = length(index);
                else
                    start(iDim) = 0;
                    count(iDim) = self.dimensions(iDim).nPoints;
                end
            end
            netcdf.putVar(self.group.id, self.id, start, count, data);
            if dimension.isMutable
                dimension.updateLength;
            end
            % netcdf.sync(self.group.id);
        end


    end
end