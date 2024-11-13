classdef NetCDFVariable < handle
    % Encapsulates variable data in a NetCDF file
    % 

    properties
        group
        id
        name
        dimensions
        attributes
        type
    end

    properties (Dependent)
        value
    end

    properties (Constant)
        % A Boolean value that indicates whether the variable is complex valued
        % - Topic: Schema keys — Variables
        GLNetCDFSchemaIsComplexKey = "isComplex";

        % A Boolean value that indicates whether this is the real part of the variable
        % - Topic: Schema keys — Variables
        GLNetCDFSchemaIsRealPartKey = "isRealPart";

        % A Boolean value that indicates whether this is the complex part of the variable
        % - Topic: Schema keys — Variables
        GLNetCDFSchemaIsImaginaryPartKey = "isImaginaryPart";

        % A Boolean value that indicates whether the variable was defined as a row vector
        % - Topic: Schema keys — Variables
        GLNetCDFSchemaIsRowVectorKey = "isRowVector";
    end


    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = NetCDFVariable(group,options)
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
            % aVariable = NetCDFVariable(group,name=aName,attributes=someAttributes,type='NC_DOUBLE');
            % ```
            %
            arguments
                group (1,1) NetCDFGroup
                options.name string
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
                self (1,1) NetCDFVariable
                group (1,1) NetCDFGroup
                id (1,1) double
            end
            self.group = group;
            self.id = id;
            self.attributes = containers.Map;

            [self.name, xtype, dimids, numatts] = netcdf.inqVar(group.id,id);
            self.type = NetCDFVariable.typeStringForTypeID(xtype);
            
            for iAtt=0:(numatts-1)
                gattname = netcdf.inqAttName(group.id,id,iAtt);
                self.attributes(gattname) = netcdf.getAtt(group.id,id,gattname);
            end
            self.dimensions = self.group.dimensionWithID(dimids);
        end

        function initVariable(self,group,name,dimensions,ncType,attributes)
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
            % ncfile.initVariable('fluid-tracer', {'x','y','t'},'NC_DOUBLE',containers.Map({'isTracer'},{'1'}));
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
                self (1,1) NetCDFVariable {mustBeNonempty}
                group (1,1) NetCDFGroup {mustBeNonempty}
                name char {mustBeNonempty}
                dimensions (:,1) NetCDFDimension
                % dimNames (:,1) cell
                ncType string {mustBeNonempty}
                attributes containers.Map = containers.Map(KeyType='double',ValueType='any')         
            end
            self.group = group;
            self.dimensions = dimensions; %group.dimensionWithName(dimNames);
            self.id = netcdf.defVar(group.id, name, ncType, [self.dimensions.id]);
            self.name = name;
            self.type = ncType;

            if ~isempty(attributes)
                keyNames = keys(attributes);
                for iKey = 1:length(keyNames)
                    netcdf.putAtt(group.id,self.id, keyNames{iKey}, attributes(keyNames{iKey}));
                end
            else
                attributes = containers.Map();
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
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Read/write data
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function value = get.value(self)
            arguments
                self NetCDFVariable
            end
            value = netcdf.getVar(self.group.id,self.id);
            if isvector(value) && isKey(self.attributes, self.GLNetCDFSchemaIsRowVectorKey)
                if self.attributes(self.GLNetCDFSchemaIsRowVectorKey) == 1
                    value=reshape(value,1,[]);
                else
                    value=reshape(value,[],1);
                end
            end
        end

        function set.value(self,data)
            arguments
                self NetCDFVariable {mustBeNonempty}
                data {mustBeNonempty}
            end
            for iDim=1:length(self.dimensions)
                if (length(self.dimensions)==1 && length(data) ~= self.dimensions(iDim).nPoints) || (length(self.dimensions)~=1 && size(data,iDim) ~= self.dimensions(iDim).nPoints)
                    error('Incorrect dimension size: dimension %d of the data of %s is length %d, but the dimension %s has length %d.',iDim,self.name,size(data,iDim),self.dimensions(iDim).name,self.dimensions(iDim).nPoints);
                end
            end
            netcdf.putVar(self.group.id, self.id, data);
            if isvector(data)
                self.addAttribute(self.GLNetCDFSchemaIsRowVectorKey, uint8(isrow(data)));
            end
        end

        function value = valueAlongDimensionAtIndex(self,dimensionName,index)
            arguments
                self NetCDFVariable {mustBeNonempty} 
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
            % variable.concatenateValueAlongDimensionAtIndex(data,dim,outputIndex);
            % ```
            %
            % - Topic: Working with variables
            % - Declaration: concatenateValueAlongDimensionAtIndex(data,dimension,index)
            % - Parameter data: variable data
            % - Parameter dimension: the variable dimension along which to concatenate
            % - Parameter index: index at which to write data
            arguments
                self NetCDFVariable {mustBeNonempty} 
                data 
                dimensionName char {mustBeNonempty}
                index {mustBeNonnegative}
            end
            dimension = self.group.dimensionWithName(dimensionName);
            if dimension.isMutable ~= 1
                error('Cannot concatenate along a dimension that is not mutable');
            end

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
            netcdf.putVar(self.group.id, self.id, start, count, data);
        end


    end

    methods (Static)
        function val = netCDFTypeForData(data)
            keys = {'double','single','int64','uint64','int32','uint32','int16','uint16','int8','uint8','char','string','logical'};
            values = {'NC_DOUBLE','NC_FLOAT','NC_INT64','NC_UINT64','NC_INT','NC_UINT','NC_SHORT','NC_USHORT','NC_BYTE','NC_UBYTE','NC_CHAR','NC_CHAR','NC_BYTE'};
            map = containers.Map(keys, values);
            if ~isKey(map,class(data))
                error('unknown data type');
            end
            val = map(class(data));
        end

        function val = netCDF3TypeForData(data)
            keys = {'double','single','int64','uint64','int32','uint32','int16','uint16','int8','uint8','char','string','logical'};
            values = {'NC_DOUBLE','NC_FLOAT','NC_INT64','NC_INT64','NC_INT','NC_INT','NC_SHORT','NC_SHORT','NC_BYTE','NC_BYTE','NC_CHAR','NC_CHAR','NC_BYTE'};
            map = containers.Map(keys, values);
            if ~isKey(map,class(data))
                error('unknown data type');
            end
            val = map(class(data));
        end

        function val = typeStringForTypeID(type)
            types = {'NC_DOUBLE','NC_FLOAT','NC_INT64','NC_UINT64','NC_INT','NC_UINT','NC_SHORT','NC_USHORT','NC_BYTE','NC_UBYTE','NC_CHAR','NC_CHAR'};
            val = nan;
            for i=1:length(types)
                if netcdf.getConstant(types{i}) == type
                    val = types{i};
                    return
                end
            end
        end
    end
end