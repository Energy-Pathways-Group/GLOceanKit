classdef NetCDFFile < handle
    %NetCDFFile Class for reading and writing to NetCDF files
    %   Simplifies and improves the consistency of reading and writing to
    %   NetCDF files compared to the built-in Matlab options.

    properties
        path    % path the NetCDF file
        ncid    % NetCDF file handle

        dimensions  % array of NetCDFDimension objects
        variables   % array of NetCDFVariable objects
        attributes  % key-value Map of global attributes

        dimensionWithName   % key-value Map to retrieve a NetCDFDimension object by name
        variableWithName    % key-value Map to retrieve a NetCDFVariable object by name

        complexVariables            % array of NetCDFComplexVariable objects
        complexVariableWithName     % key-value Map to retrieve a NetCDFComplexVariable object by name
    end

    properties (Constant)
        GLNetCDFSchemaVersionKey = "GLNetCDFSchemaVersion";

        % attributes for coordinate variables (dimensions)
        GLNetCDFSchemaIsCoordinateVariableKey = "isCoordinateVariable";
        GLNetCDFSchemaIsPeridiocKey = "isPeriodic";
        GLNetCDFSchemaMutableKey = "isMutable";
        GLNetCDFSchemaDomainMinimumKey = "domainMin";
        GLNetCDFSchemaBasisFunctionKey = "basisFunction";
        GLNetCDFSchemaDomainLengthKey = "domainLength";
        GLNetCDFSchemaIsEvenlySampledKey = "isEvenlySampled";
        GLNetCDFSchemaSampleIntervalKey = "sampleInterval";
        GLNetCDFSchemaGridTypeKey = "gridType";
        GLNetCDFSchemaIsFrequencyDomainKey = "isFrequencyDomain";

        % attributes for variables and dimensions
        GLNetCDFSchemaUnitsKey = "units";

        % attributes for variables
        GLNetCDFSchemaIsComplexKey = "isComplex";
        GLNetCDFSchemaProperNameKey = "properName";
        GLNetCDFSchemaIsRealPartKey = "isRealPart";
        GLNetCDFSchemaIsImaginaryPartKey = "isImaginaryPart";
        GLNetCDFSchemaIsRowVectorKey = "isRowVector";
        
        GLNetCDFSchemaUniqueVariableIDKey = "uniqueVariableID";
    end

    methods
        function self = NetCDFFile(path,overwriteExisting)
            % NetCDFFile initialize an from existing or create new file
            %
            % Calling,
            %   ncfile = NetCDFFile(path)
            % will load an existing file (if one exists) or create a new
            % file (if none exists).
            %
            %   ncfile = NetCDFFile(path,'OVERWRITE_EXISTING')
            % will delete any existing file and create a new file.

            self.path = path;
            shouldOverwrite = 0;
            self.dimensions = NetCDFDimension.empty(0,0);
            self.dimensionWithName = containers.Map;
            self.variables = NetCDFVariable.empty(0,0);
            self.variableWithName = containers.Map;
            self.attributes = containers.Map;
            self.complexVariables = NetCDFComplexVariable.empty(0,0);
            self.complexVariableWithName = containers.Map;

            if nargin == 2 && strcmp(overwriteExisting,'OVERWRITE_EXISTING')
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

        function complexVariable = initComplexVariable(self,name,dimNames,ncType,properties)
            if isKey(self.complexVariableWithName,name)
                error('A variable with that name already exists.');
            end
            if nargin < 5 || isempty(properties)
                properties = containers.Map;
            end
            properties(self.GLNetCDFSchemaIsComplexKey) = 1;

            properties(self.GLNetCDFSchemaIsRealPartKey) = 1;
            properties(self.GLNetCDFSchemaIsImaginaryPartKey) = 0;
            realVar = self.initVariable(strcat(name,"_realp"),dimNames,ncType,properties);

            properties(self.GLNetCDFSchemaIsRealPartKey) = 0;
            properties(self.GLNetCDFSchemaIsImaginaryPartKey) = 1;
            imagVar = self.initVariable(strcat(name,"_imagp"),dimNames,ncType,properties);
            
            complexVariable = NetCDFComplexVariable(name,realVar,imagVar);
            self.complexVariables(end+1) = complexVariable;
            self.complexVariableWithName(name) = complexVariable;
        end

        function variable = initVariable(self,name,dimNames,ncType,properties)
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

        function setVariable(self,data,name)
            if isKey(self.complexVariableWithName,name)
                complexVariable = self.complexVariableWithName(name);
                variable = complexVariable.realVar;
                for iDim=1:length(variable.dimensions)
                    if (isvector(data) && length(data) ~= variable.dimensions(iDim).nPoints) || (~isvector(data) && size(data,iDim) ~= variable.dimensions(iDim).nPoints)
                        error('Incorrect dimension size: dimension %d of the data is length %d, but the dimension %s has length %d.',iDim,size(data,iDim),variable.dimensions(iDim).name,variable.dimensions(iDim).nPoints);
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
                        error('Incorrect dimension size: dimension %d of the data is length %d, but the dimension %s has length %d.',iDim,size(data,iDim),variable.dimensions(iDim).name,variable.dimensions(iDim).nPoints);
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
            if nargin < 6 || isempty(ncType)
                ncType = self.netCDFTypeForData(data);
            end

            if isreal(data)
                variable = self.initVariable(name,dimNames,ncType,properties);
            else
                variable = self.initComplexVariable(name,dimNames,ncType,properties);
            end

            self.setVariable(data,name);
        end

        function concatenateVariableAlongDimension(self,name,data,dimName,index)
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
            dims = NetCDFDimension.empty(length(dimids),0);
            for iDim=1:length(dimids)
                dims(iDim) = self.dimensions(dimids(iDim)+1);
            end
        end

        function val = netCDFTypeForData(self,data)
            keys = {'double','single','int64','uint64','int32','uint32','int16','uint16','int8','uint8','char','string'};
            values = {'NC_DOUBLE','NC_FLOAT','NC_INT64','NC_UINT64','NC_INT','NC_UINT','NC_SHORT','NC_USHORT','NC_BYTE','NC_UBYTE','NC_CHAR','NC_CHAR'};
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
            netcdf.sync(self.ncid);
        end

        function self = open(self)
            self.ncid = netcdf.open(self.path, bitor(netcdf.getConstant('SHARE'),netcdf.getConstant('WRITE')));
        end

        function self = close(self)
            netcdf.close(self.ncid);
            self.ncid = [];
        end
    end
end