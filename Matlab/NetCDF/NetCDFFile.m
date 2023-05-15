classdef NetCDFFile < handle
    % A class for reading and writing to NetCDF files
    %
    % NetCDF files are a standard file format for reading and writing data.
    % This class is designed to simplify the task of adding new dimensions,
    % variables, and attributes to a NetCDF file compared to using the
    % built-in `ncread` and `ncwrite` functions.
    %
    % ```matlab
    % ncfile = NetCDFFile('myfile.nc')
    %
    % % create two new dimensions and add them to the file
    % x = linspace(0,10,11);
    % y = linspace(-10,0,11);
    % ncfile.addDimension('x',x);
    % ncfile.addDimension('y',y);
    %
    % % Create new multi-dimensional variables, and add those to the file
    % [X,Y] = ncgrid(x,y);
    % ncfile.addVariable(X,{'x','y'});
    % ncfile.addVariable(Y,{'x','y'});
    % ```
    %
    % - Topic: Initializing
    % - Topic: Accessing file properties
    % - Topic: Working with dimensions
    % - Topic: Working with variables
    % - Topic: Working with global attributes
    % - Topic: Schema keys
    % - Topic: Schema keys — Dimensions
    % - Topic: Schema keys — Variables
    %
    % - Declaration: classdef NetCDFFile < handle
    properties
        % file path the NetCDF file
        % - Topic: Accessing file properties
        path

        % file handle
        % - Topic: Accessing file properties
        ncid

        % array of NetCDFDimension objects
        %
        % An array of NetCDFDimension objects for each coordinate dimension
        % defined in the NetCDF file. The dimensions order in the array
        % should reflect the underlying dimensionID defined in the NetCDF
        % file.
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

        % key-value Map to retrieve a NetCDFDimension object by name
        %
        % Usage
        % ```matlab
        % xDim = ncfile.dimensionWithName('x');
        % ```
        %
        % - Topic: Working with dimensions
        dimensionWithName   

        % key-value Map to retrieve a NetCDFVariable object by name
        % - Topic: Working with variables
        variableWithName

        % array of NetCDFComplexVariable objects
        % - Topic: Working with variables
        complexVariables            

        % key-value Map to retrieve a NetCDFComplexVariable object by name
        % - Topic: Working with variables
        complexVariableWithName     
    end

    properties (Constant)
        % - Topic: Schema keys
        GLNetCDFSchemaVersionKey = "GLNetCDFSchemaVersion";

        % A Boolean value that indicates whether the dimension is associated with a coordinate variable
        % - Topic: Schema keys — Dimensions
        GLNetCDFSchemaIsCoordinateVariableKey = "isCoordinateVariable";

        % A Boolean value that indicates whether the dimension is periodic
        % - Topic: Schema keys — Dimensions
        GLNetCDFSchemaIsPeridiocKey = "isPeriodic";

        % A Boolean value that indicates whether the dimension is mutable
        % - Topic: Schema keys — Dimensions
        GLNetCDFSchemaMutableKey = "isMutable";

        % The minimum value of the domain
        % - Topic: Schema keys — Dimensions
        GLNetCDFSchemaDomainMinimumKey = "domainMin";

        % What basis function describe this dimension
        % - Topic: Schema keys — Dimensions
        GLNetCDFSchemaBasisFunctionKey = "basisFunction";

        % The length of the domain
        % - Topic: Schema keys — Dimensions
        GLNetCDFSchemaDomainLengthKey = "domainLength";

        % A Boolean value that indicates whether the dimension has even sampling
        % - Topic: Schema keys — Dimensions
        GLNetCDFSchemaIsEvenlySampledKey = "isEvenlySampled";

        % sample interval of the domain, if it is evenly sampled
        % - Topic: Schema keys — Dimensions
        GLNetCDFSchemaSampleIntervalKey = "sampleInterval";

        % type of grid
        % - Topic: Schema keys — Dimensions
        GLNetCDFSchemaGridTypeKey = "gridType";

        % A Boolean value that indicates whether the dimension is considered in the frequency (spectral) domain
        % - Topic: Schema keys — Dimensions
        GLNetCDFSchemaIsFrequencyDomainKey = "isFrequencyDomain";

        % Units of the variable or dimension
        % - Topic: Schema keys
        GLNetCDFSchemaUnitsKey = "units";

        % A Boolean value that indicates whether the variable is complex valued
        % - Topic: Schema keys — Variables
        GLNetCDFSchemaIsComplexKey = "isComplex";

        % Human readable name of the variable
        % - Topic: Schema keys — Variables
        GLNetCDFSchemaProperNameKey = "properName";

        % A Boolean value that indicates whether this is the real part of the variable
        % - Topic: Schema keys — Variables
        GLNetCDFSchemaIsRealPartKey = "isRealPart";

        % A Boolean value that indicates whether this is the complex part of the variable
        % - Topic: Schema keys — Variables
        GLNetCDFSchemaIsImaginaryPartKey = "isImaginaryPart";

        % A Boolean value that indicates whether the variable was defined as a row vector
        % - Topic: Schema keys — Variables
        GLNetCDFSchemaIsRowVectorKey = "isRowVector";

        % A custom unique variable ID
        % - Topic: Schema keys — Variables
        GLNetCDFSchemaUniqueVariableIDKey = "uniqueVariableID";
    end

    methods
        function self = NetCDFFile(path,options)
            % NetCDFFile initialize an from existing or create new file
            %
            % Calling,
            %   ncfile = NetCDFFile(path)
            % will load an existing file (if one exists) or create a new
            % file (if none exists).
            %
            %   ncfile = NetCDFFile(path,shouldOverwriteExisting=1)
            % will delete any existing file and create a new file.
            %
            % - Topic: Initialization
            % - Declaration: ncfile = etCDFFile(path,options)
            % - Parameter path: path to write file
            % - Parameter shouldOverwriteExisting: (optional) boolean indicating whether or not to overwrite an existing file at the path. Default 0.
            % - Returns: a new NetCDFFile instance
            arguments
                path char {mustBeNonempty}
                options.shouldOverwriteExisting double {mustBeMember(options.shouldOverwriteExisting,[0 1])} = 0 
            end
            self.path = path;
            shouldOverwrite = 0;
            self.dimensions = NetCDFDimension.empty(0,0);
            self.dimensionWithName = containers.Map;
            self.variables = NetCDFVariable.empty(0,0);
            self.variableWithName = containers.Map;
            self.attributes = containers.Map;
            self.complexVariables = NetCDFComplexVariable.empty(0,0);
            self.complexVariableWithName = containers.Map;

            if options.shouldOverwriteExisting == 1
                shouldOverwrite = 1;
            end
            if isfile(self.path)
                if shouldOverwrite == 1
                    delete(self.path);
                    self.CreateNewFile();
                else
                    self.InitializeFromExistingFile();
                end
            else
                self.CreateNewFile();
            end
        end

        function CreateNewFile(self)
            self.ncid = netcdf.create(self.path, bitor(netcdf.getConstant('SHARE'),netcdf.getConstant('WRITE')));
            netcdf.endDef(self.ncid);
        end

        function InitializeFromExistingFile(self)
            self.open();

            [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(self.ncid);
            
            % Fetch all dimensions and convert to NetCDFDimension objects
            for iDim=0:(ndims-1)
                [dimname,dimlen] = netcdf.inqDim(self.ncid,iDim);
                self.dimensions(iDim+1) = NetCDFDimension(dimname,dimlen,iDim);
                self.dimensionWithName(dimname) = self.dimensions(iDim+1);
                if iDim == unlimdimid
                    self.dimensions(iDim+1).isMutable = 1;
                end
            end
            self.dimensions = self.dimensions.';

            % Fetch all variables and convert to NetCDFVariable objects
            for iVar=0:(nvars-1)
                [varname, xtype, dimids, numatts] = netcdf.inqVar(self.ncid,iVar);
                variableAttributes = containers.Map;
                for iAtt=0:(numatts-1)
                    gattname = netcdf.inqAttName(self.ncid,iVar,iAtt);
                    variableAttributes(gattname) = netcdf.getAtt(self.ncid,iVar,gattname);
                end
                self.variables(iVar+1) = NetCDFVariable(varname,self.dimensionsForDimIDs(dimids),variableAttributes,self.typeStringForTypeID(xtype),iVar);
                self.variableWithName(varname) = self.variables(iVar+1);
            end
            self.variables = self.variables.';
 
            % Grab all the global attributes
            for iAtt=0:(ngatts-1)
                gattname = netcdf.inqAttName(self.ncid,netcdf.getConstant('NC_GLOBAL'),iAtt);
                self.attributes(gattname) = netcdf.getAtt(self.ncid,netcdf.getConstant('NC_GLOBAL'),gattname);
            end

            % Now check to see if any variables form a complex variable
            for iVar=1:length(self.variables)
                if isKey(self.variables(iVar).attributes,{self.GLNetCDFSchemaIsComplexKey,self.GLNetCDFSchemaIsRealPartKey})
                    if self.variables(iVar).attributes(self.GLNetCDFSchemaIsComplexKey) == 1 && self.variables(iVar).attributes(self.GLNetCDFSchemaIsRealPartKey) == 1
                        complexName = extractBefore(self.variables(iVar).name,"_realp");
                        imagName = strcat(complexName,"_imagp");
                        if isKey(self.variableWithName,imagName)
                            complexVariable = NetCDFComplexVariable(complexName,self.variables(iVar),self.variableWithName(imagName));
                            self.complexVariables(end+1) = complexVariable;
                            self.complexVariableWithName(complexName) = complexVariable;
                        end
                    end
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Attributes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        function addAttribute(self,name,data)
            % add a global attribute to the file
            %
            % - Topic: Working with global attributes
            % - Declaration: addAttribute(name,data)
            % - Parameter name: string of the attribute name
            % - Parameter data: value
            netcdf.reDef(self.ncid);
            netcdf.putAtt(self.ncid,netcdf.getConstant('NC_GLOBAL'), name, data);
            netcdf.endDef(self.ncid);
            self.attributes(name) = data;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Dimensions
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [dimension,variable] = addDimension(self,name,data,properties,dimLength)
            % Adds a both a new dimension and its associated coordinate variable to the NetCDF file.
            %
            % Usage
            % ```matlab
            % x = linspace(0,10,11);
            % ncfile.addDimension('x',x,[]);
            % ```
            %
            % - Topic: Working with dimensions
            % - Declaration: [dimension,variable] = addDimension(name,data,properties,dimLength)
            % - Parameter name: string with the name of the dimension
            % - Parameter data: array of values along that dimension, or empty
            % - Parameter properties: containers.Map containing any key-value pairs to be associated with the dimension.
            % - Parameter dimLength: (optional) length of the dimension
            % - Returns dimension: a NetCDFDimension object with the newly create dimension
            % - Returns variable: a NetCDFVariable object with the associated coordinate variable
            if nargin < 5
                dimLength = 0;
            end
            if isKey(self.dimensionWithName,name) || isKey(self.variableWithName,name)
                error('A dimension with that name already exists.');
            end
            if (dimLength == 0 && isempty(data)) || isempty(name)
                error('You must specify a name and data');
            end
            if isinf(dimLength)
                n = netcdf.getConstant('NC_UNLIMITED');
            elseif isempty(data)
                n = dimLength;
            else
                n = length(data);
            end
            
            ncType = self.netCDFTypeForData(data);
            netcdf.reDef(self.ncid);
            dimID = netcdf.defDim(self.ncid, name, n);
            varID = netcdf.defVar(self.ncid, name, ncType, dimID);
            if ~isempty(properties)
                keyNames = keys(properties);
                for iKey = 1:length(keyNames)
                    netcdf.putAtt(self.ncid,varID, keyNames{iKey}, properties(keyNames{iKey}));
                end
            end
            netcdf.endDef(self.ncid);
            if ~isempty(data)
                netcdf.putVar(self.ncid, varID, data);
            end
            
            dimension = NetCDFDimension(name,n,dimID);
            variable = NetCDFVariable(name,dimension,properties,ncType,varID);

            % Can this fail? I don't think it should.
            self.dimensions(dimID+1) = dimension;
            self.dimensionWithName(name) = dimension;
            self.variables(varID+1) = variable;
            self.variableWithName(name) = variable;
        end

        function [dimension,variable] = addMutableDimension(self,name,properties)
            [dimension,variable] = self.addDimension(name,[],properties,inf);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Variables
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %addVariable(self,name,data,dimNames,properties,ncType)
        function complexVariable = initComplexVariable(self,name,dimNames,properties,ncType)
            % initialize a complex-valued variable
            %
            % NetCDF does not directly work with complex variables, so this
            % method manages the hassle of working with the real and
            % imaginary parts separately.
            %
            % - Topic: Working with variables
            % - Declaration: complexVariable = initComplexVariable(name,dimNames,properties,ncType)
            % - Parameter name: name of the variable (a string)
            % - Parameter dimNames: cell array containing the dimension names
            % - Parameter properties: (optional) `containers.Map`
            % - Parameter properties: ncType
            % - Returns complexVariable: a NetCDFComplexVariable object
            if isKey(self.complexVariableWithName,name)
                error('A variable with the name %s already exists.',name);
            end
            if isempty(properties)
                properties = containers.Map;
            end
            properties(self.GLNetCDFSchemaIsComplexKey) = 1;

            properties(self.GLNetCDFSchemaIsRealPartKey) = 1;
            properties(self.GLNetCDFSchemaIsImaginaryPartKey) = 0;
            realVar = self.initVariable(strcat(name,"_realp"),dimNames,properties,ncType);

            properties(self.GLNetCDFSchemaIsRealPartKey) = 0;
            properties(self.GLNetCDFSchemaIsImaginaryPartKey) = 1;
            imagVar = self.initVariable(strcat(name,"_imagp"),dimNames,properties,ncType);
            
            complexVariable = NetCDFComplexVariable(name,realVar,imagVar);
            self.complexVariables(end+1) = complexVariable;
            self.complexVariableWithName(name) = complexVariable;
        end

        function variable = initVariable(self,name,dimNames,properties,ncType)
            % initialize a real-valued variable
            %
            %
            % - Topic: Working with variables
            if isKey(self.variableWithName,name)
                error('A variable with that name already exists.');
            end
            if ischar(dimNames)
                dimNames = {dimNames};
            end
            dimIDs = nan(length(dimNames),1);
            for iDim=1:length(dimNames)
                if ~isKey(self.dimensionWithName,dimNames{iDim})
                    error('Unable to find a dimension with the name %s',dimNames{iDim});
                end
                dimIDs(iDim) = self.dimensionWithName(dimNames{iDim}).dimID;
            end

            netcdf.reDef(self.ncid);
            varID = netcdf.defVar(self.ncid, name, ncType, dimIDs);
            if ~isempty(properties)
                keyNames = keys(properties);
                for iKey = 1:length(keyNames)
                    netcdf.putAtt(self.ncid,varID, keyNames{iKey}, properties(keyNames{iKey}));
                end
            end
            netcdf.endDef(self.ncid);

            variable = NetCDFVariable(name,self.dimensionsForDimIDs(dimIDs),properties,ncType,varID);
            self.variables(varID+1) = variable;
            self.variableWithName(name) = variable;
        end

        function setVariable(self,name,data)
            % add data for a variable with a given name
            %
            % - Topic: Working with variables
            if ~isvector(data) && length(data) == numel(data)
                % rare case that something is, e.g., size(data) = [1 1 N]
                data=squeeze(data);
            end
            if isKey(self.complexVariableWithName,name)
                complexVariable = self.complexVariableWithName(name);
                variable = complexVariable.realVar;
                for iDim=1:length(variable.dimensions)
                    if (isvector(data) && length(data) ~= variable.dimensions(iDim).nPoints) || (~isvector(data) && size(data,iDim) ~= variable.dimensions(iDim).nPoints)
                        error('Incorrect dimension size: dimension %d of the data of %s is length %d, but the dimension %s has length %d.',iDim,name,size(data,iDim),variable.dimensions(iDim).name,variable.dimensions(iDim).nPoints);
                    end
                end
                netcdf.putVar(self.ncid, complexVariable.realVar.varID, real(data));
                netcdf.putVar(self.ncid, complexVariable.imagVar.varID, imag(data));
                if isvector(data)
                    complexVariable.realVar.attributes(self.GLNetCDFSchemaIsRowVectorKey) = isrow(data);
                    complexVariable.imagVar.attributes(self.GLNetCDFSchemaIsRowVectorKey) = isrow(data);
                    netcdf.reDef(self.ncid);
                    netcdf.putAtt(self.ncid,complexVariable.realVar.varID, self.GLNetCDFSchemaIsRowVectorKey, uint8(isrow(data)));
                    netcdf.putAtt(self.ncid,complexVariable.imagVar.varID, self.GLNetCDFSchemaIsRowVectorKey, uint8(isrow(data)));
                    netcdf.endDef(self.ncid);
                end
            else
                variable = self.variableWithName(name);
                for iDim=1:length(variable.dimensions)
                    if (isvector(data) && length(data) ~= variable.dimensions(iDim).nPoints) || (~isvector(data) && size(data,iDim) ~= variable.dimensions(iDim).nPoints)
                        error('Incorrect dimension size: dimension %d of the data of %s is length %d, but the dimension %s has length %d.',iDim,name,size(data,iDim),variable.dimensions(iDim).name,variable.dimensions(iDim).nPoints);
                    end
                end
                netcdf.putVar(self.ncid, variable.varID, data);
                if isvector(data)
                    variable.attributes(self.GLNetCDFSchemaIsRowVectorKey) = isrow(data);
                    netcdf.reDef(self.ncid);
                    netcdf.putAtt(self.ncid,variable.varID, self.GLNetCDFSchemaIsRowVectorKey, uint8(isrow(data)));
                    netcdf.endDef(self.ncid);
                end
            end

        end

        function varargout = readVariables(self,varargin)
            % read data from variables
            %
            % - Topic: Working with variables
            varargout = cell(size(varargin));
            for iArg=1:length(varargin)
                if isKey(self.complexVariableWithName,varargin{iArg})
                    varargout{iArg} = complex(netcdf.getVar(self.ncid,self.complexVariableWithName(varargin{iArg}).realVar.varID),netcdf.getVar(self.ncid,self.complexVariableWithName(varargin{iArg}).imagVar.varID));
                    if isvector(varargout{iArg}) && isKey(self.complexVariableWithName(varargin{iArg}).realVar.attributes, self.GLNetCDFSchemaIsRowVectorKey)
                        if self.complexVariableWithName(varargin{iArg}).realVar.attributes(self.GLNetCDFSchemaIsRowVectorKey) == 1
                            varargout{iArg}=reshape(varargout{iArg},1,[]);
                        else
                            varargout{iArg}=reshape(varargout{iArg},[],1);
                        end
                    end
                elseif isKey(self.variableWithName,varargin{iArg})
                    varargout{iArg} = netcdf.getVar(self.ncid,self.variableWithName(varargin{iArg}).varID);
                    if isvector(varargout{iArg}) && isKey(self.variableWithName(varargin{iArg}).attributes, self.GLNetCDFSchemaIsRowVectorKey)
                        if self.variableWithName(varargin{iArg}).attributes(self.GLNetCDFSchemaIsRowVectorKey) == 1
                            varargout{iArg}=reshape(varargout{iArg},1,[]);
                        else
                            varargout{iArg}=reshape(varargout{iArg},[],1);
                        end
                    end
                else
                    error('Unable to find a variable with the name %s',varargin{iArg});
                end
            end
        end

        function varargout = readVariablesAtIndexAlongDimension(self,dimName,index,varargin)
            % read data from variables from a particular index
            %
            % - Topic: Working with variables
            varargout = cell(size(varargin));
            for iArg=1:length(varargin)
                if isKey(self.complexVariableWithName,varargin{iArg})
                    complexVariable = self.complexVariableWithName(varargin{iArg});
                    variable = complexVariable.realVar;
                    start = zeros(1,length(variable.dimensions));
                    count = zeros(1,length(variable.dimensions));
                    for iDim=1:length(variable.dimensions)
                        if strcmp(variable.dimensions(iDim).name,dimName)
                            start(iDim) = index-1;
                            count(iDim) = 1;
                        else
                            start(iDim) = 0;
                            count(iDim) = variable.dimensions(iDim).nPoints;
                        end
                    end
                varargout{iArg} = complex(netcdf.getVar(self.ncid, complexVariable.realVar.varID, start, count),netcdf.getVar(self.ncid, complexVariable.imagVar.varID, start, count));
                elseif isKey(self.variableWithName,varargin{iArg})
                    variable = self.variableWithName(varargin{iArg});
                    start = zeros(1,length(variable.dimensions));
                    count = zeros(1,length(variable.dimensions));
                    for iDim=1:length(variable.dimensions)
                        if strcmp(variable.dimensions(iDim).name,dimName)
                            start(iDim) = index-1;
                            count(iDim) = 1;
                        else
                            start(iDim) = 0;
                            count(iDim) = variable.dimensions(iDim).nPoints;
                        end
                    end
                    varargout{iArg} = netcdf.getVar(self.ncid, variable.varID, start, count);
                else
                    error('Unable to find a variable with the name %s',name);
                end
            end
        end

        function variable = addVariable(self,name,data,dimNames,properties,ncType)
            % add a new variable to the file
            %
            % - Topic: Working with variables
            if nargin < 6 || isempty(ncType)
                ncType = self.netCDFTypeForData(data);
            end
            if nargin < 5
                properties = [];
            end

            if isreal(data)
                variable = self.initVariable(name,dimNames,properties,ncType);
            else
                variable = self.initComplexVariable(name,dimNames,properties,ncType);
            end

            self.setVariable(name,data);
        end

        function concatenateVariableAlongDimension(self,name,data,dimName,index)
            % append new data to an existing variable
            %
            % - Topic: Working with variables
            if isKey(self.complexVariableWithName,name)
                complexVariable = self.complexVariableWithName(name);
                variable = complexVariable.realVar;
                start = zeros(1,length(variable.dimensions));
                count = zeros(1,length(variable.dimensions));
                for iDim=1:length(variable.dimensions)
                    if strcmp(variable.dimensions(iDim).name,dimName)
                        start(iDim) = index-1;
                        count(iDim) = 1;
                    else
                        start(iDim) = 0;
                        count(iDim) = variable.dimensions(iDim).nPoints;
                    end
                end
                netcdf.putVar(self.ncid, complexVariable.realVar.varID, start, count, real(data));
                netcdf.putVar(self.ncid, complexVariable.imagVar.varID, start, count, imag(data));
            elseif isKey(self.variableWithName,name)
                variable = self.variableWithName(name);
                start = zeros(1,length(variable.dimensions));
                count = zeros(1,length(variable.dimensions));
                for iDim=1:length(variable.dimensions)
                    if strcmp(variable.dimensions(iDim).name,dimName)
                        start(iDim) = index-1;
                        count(iDim) = 1;
                    else
                        start(iDim) = 0;
                        count(iDim) = variable.dimensions(iDim).nPoints;
                    end
                end
                netcdf.putVar(self.ncid, variable.varID, start, count, data);
            else
                error('Unable to find a variable with the name %s',name);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Utilities
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function dims = dimensionsForDimIDs(self,dimids)
            % return the dimension IDs given the dimension names
            %
            % - Topic: Working with dimensions
            dims = NetCDFDimension.empty(length(dimids),0);
            for iDim=1:length(dimids)
                dims(iDim) = self.dimensions(dimids(iDim)+1);
            end
        end

        function val = netCDFTypeForData(self,data)
            keys = {'double','single','int64','uint64','int32','uint32','int16','uint16','int8','uint8','char','string','logical'};
            values = {'NC_DOUBLE','NC_FLOAT','NC_INT64','NC_UINT64','NC_INT','NC_UINT','NC_SHORT','NC_USHORT','NC_BYTE','NC_UBYTE','NC_CHAR','NC_CHAR','NC_BYTE'};
            map = containers.Map(keys, values);
            if ~isKey(map,class(data))
                error('unknown data type');
            end
            val = map(class(data));
        end

        function val = typeStringForTypeID(self,type)
            types = {'NC_DOUBLE','NC_FLOAT','NC_INT64','NC_UINT64','NC_INT','NC_UINT','NC_SHORT','NC_USHORT','NC_BYTE','NC_UBYTE','NC_CHAR','NC_CHAR'};
            val = nan;
            for i=1:length(types)
                if netcdf.getConstant(types{i}) == type
                    val = types{i};
                    return
                end
            end
        end

        function self = sync(self)
            % - Topic: Accessing file properties
            netcdf.sync(self.ncid);
        end

        function self = open(self)
            % - Topic: Accessing file properties
            self.ncid = netcdf.open(self.path, bitor(netcdf.getConstant('SHARE'),netcdf.getConstant('WRITE')));
        end

        function self = close(self)
            % - Topic: Accessing file properties
            netcdf.close(self.ncid);
            self.ncid = [];
        end
    end
end