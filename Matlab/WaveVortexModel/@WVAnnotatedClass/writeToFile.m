function ncfile = writeToFile(self,path,variables,options)
    % Output the `WVTransform` to file.
    %
    % Writes the WVTransform instance to file, with enough information to
    % re-initialize. Pass additional variables to the variable list that
    % should also be written to file.
    %
    % Subclasses should add any necessary properties or variables to the
    % variable list before calling this superclass method.
    %
    % - Topic: Write to file
    % - Declaration: ncfile = writeToFile(netcdfFile,variables,options)
    % - Parameter path: path to write file
    % - Parameter variables: strings of variable names.
    % - Parameter shouldOverwriteExisting: (optional) boolean indicating whether or not to overwrite an existing file at the path. Default 0. 
    % - Parameter shouldAddDefaultVariables: (optional) boolean indicating whether or not add default variables `A0`,`Ap`,`Am`,`t`. Default 1.
    arguments (Input)
        self WVAnnotatedClass {mustBeNonempty}
        path char {mustBeNonempty}
    end
    arguments (Input,Repeating)
        variables char
    end
    arguments (Input)
        options.shouldOverwriteExisting double {mustBeMember(options.shouldOverwriteExisting,[0 1])} = 0
        options.shouldAddDefaultDimensions double {mustBeMember(options.shouldAddDefaultDimensions,[0 1])} = 1 
        options.shouldAddDefaultVariables double {mustBeMember(options.shouldAddDefaultVariables,[0 1])} = 1 
        options.dimensions {mustBeA(options.dimensions,"cell")} = {}
        options.attributes = configureDictionary("string","string")
    end
    arguments (Output)
        ncfile NetCDFFile
    end

    if options.shouldOverwriteExisting == 1
        if isfile(path)
            delete(path);
        end
    else
        if isfile(path)
            error('A file already exists with that name.')
        end
    end
    ncfile = NetCDFFile(path);

    if options.shouldAddDefaultDimensions == 1
        dims = union(options.dimensions,self.requiredDimensions);
    else
        dims = options.dimensions;
    end
    for iDim=1:length(dims)
        dimAnnotation = self.dimensionAnnotationWithName(dims{iDim});
        dimAnnotation.attributes('units') = dimAnnotation.units;
        dimAnnotation.attributes('long_name') = dimAnnotation.description;
        ncfile.addDimension(dimAnnotation.name,self.(dimAnnotation.name),attributes=dimAnnotation.attributes);
    end

    attributeNames = keys(options.attributes);
    for iKey=1:length(attributeNames)
        ncfile.addAttribute(attributeNames{iKey},options.attributes(attributeNames{iKey}));
    end

    if options.shouldAddDefaultVariables == 1
        variables = union(variables,self.requiredVariables);
    end
    for iVar=1:length(variables)
        if isKey(self.propertyAnnotationNameMap,variables{iVar})
            varAnnotation = self.propertyAnnotationNameMap(variables{iVar});
        elseif isKey(self.variableAnnotationWithName,variables{iVar})
            varAnnotation = self.variableAnnotationWithName(variables{iVar});
        end
        % check for function_handle, and add group
        varAnnotation.attributes('units') = varAnnotation.units;
        varAnnotation.attributes('long_name') = varAnnotation.description;
        ncfile.addVariable(varAnnotation.name,varAnnotation.dimensions,self.(varAnnotation.name),isComplex=varAnnotation.isComplex,attributes=varAnnotation.attributes);
    end
end