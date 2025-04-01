classdef NetCDFVariable < handle
    % Encapsulates variable data in a NetCDF file
    % 
    properties (WeakHandle)
        group NetCDFGroup
    end
    
    properties
        name
        dimensions
        attributes
        type
        namePath
    end

    properties (Abstract)
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

        % A Boolean value that indicates whether the variable was defined as a logical type
        % - Topic: Schema keys — Variables
        GLNetCDFSchemaIsLogicalTypeKey = "isLogicalType";

        % A Boolean value that indicates whether the variable was defined to hold a function_handle type
        % - Topic: Schema keys — Variables
        GLNetCDFSchemaIsFunctionHandleTypeKey = "isFunctionHandleType";
    end


    methods (Abstract)
        addAttribute(self,key,value)
        value = valueAlongDimensionAtIndex(self,dimensionName,index)
        setValueAlongDimensionAtIndex(self,data,dimensionName,index)
    end

    methods
        function val = get.namePath(self)
            if isequal(self.group.groupPath,"")
                val = self.name;
            else
                val = self.group.groupPath + "/" + self.name;
            end
        end
    end

    methods (Static)
        function val = matlabTypeForNetCDFType(ncType)
            % https://www.mathworks.com/help/matlab/ref/netcdf.getvar.html
            % plus the logical type, which we assign to an unsigned byte
            values = ["double" "single" "int64" "uint64" "int32" "uint32" "int16" "uint16" "int8" "uint8" "char" "string" "logical"];
            keys = ["NC_DOUBLE" "NC_FLOAT" "NC_INT64" "NC_UINT64" "NC_INT" "NC_UINT" "NC_SHORT" "NC_USHORT" "NC_BYTE" "NC_UBYTE" "NC_CHAR" "NC_STRING" "NC_UBYTE"];
            dict = dictionary(keys,values);
            val = dict(ncType);
        end

        function val = ncTypeForMatlabType(type)
            % https://www.mathworks.com/help/matlab/ref/netcdf.getvar.html
            % plus the logical type, which we assign to an unsigned byte
            keys = ["double" "single" "int64" "uint64" "int32" "uint32" "int16" "uint16" "int8" "uint8" "char" "string" "logical"];
            values = ["NC_DOUBLE" "NC_FLOAT" "NC_INT64" "NC_UINT64" "NC_INT" "NC_UINT" "NC_SHORT" "NC_USHORT" "NC_BYTE" "NC_UBYTE" "NC_CHAR" "NC_STRING" "NC_UBYTE"];
            dict = dictionary(keys,values);
            val = dict(type);
        end

        function val = netCDFTypeForData(data)
            keys = {'double','single','int64','uint64','int32','uint32','int16','uint16','int8','uint8','char','string','logical'};
            values = {'NC_DOUBLE','NC_FLOAT','NC_INT64','NC_UINT64','NC_INT','NC_UINT','NC_SHORT','NC_USHORT','NC_BYTE','NC_UBYTE','NC_CHAR','NC_STRING','NC_BYTE'};
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