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
        self PMAnnotatedClass {mustBeNonempty}
        path char {mustBeNonempty}
    end
    arguments (Input,Repeating)
        variables char
    end
    arguments (Input)
        options.shouldOverwriteExisting double {mustBeMember(options.shouldOverwriteExisting,[0 1])} = 0
        options.shouldAddRequiredDimensions double {mustBeMember(options.shouldAddRequiredDimensions,[0 1])} = 1 
        options.shouldAddRequiredVariables double {mustBeMember(options.shouldAddRequiredVariables,[0 1])} = 1 
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

    if options.shouldAddRequiredDimensions == 1
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
    ncfile.addAttribute('PMAnnotatedClass',class(self));

    attributeNames = keys(options.attributes);
    for iKey=1:length(attributeNames)
        ncfile.addAttribute(attributeNames{iKey},options.attributes(attributeNames{iKey}));
    end

    if options.shouldAddRequiredVariables == 1
        variables = union(variables,self.requiredVariables);
    end
    for iVar=1:length(variables)
        % check for function_handle, and add group
        varAnnotation = self.variableAnnotationWithName(variables{iVar});
        if ~isempty(varAnnotation.units)
            varAnnotation.attributes('units') = varAnnotation.units;
        end
        if ~isempty(varAnnotation.description)
            varAnnotation.attributes('long_name') = varAnnotation.description;
        end
        if isa(self.(varAnnotation.name),'function_handle')
            ncfile.addFunctionHandle(varAnnotation.name,self.(varAnnotation.name),attributes=varAnnotation.attributes);
        else
            ncfile.addVariable(varAnnotation.name,varAnnotation.dimensions,self.(varAnnotation.name),isComplex=varAnnotation.isComplex,attributes=varAnnotation.attributes);
        end
    end
end