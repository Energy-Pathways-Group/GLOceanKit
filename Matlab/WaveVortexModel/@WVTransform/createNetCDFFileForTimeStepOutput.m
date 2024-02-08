function [ncfile,matFilePath] = createNetCDFFileForTimeStepOutput(wvt,path,variables,options)
    % Output the `WVTransform` to file with variable time dimension
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
    arguments
        wvt WVTransform {mustBeNonempty}
        path char {mustBeNonempty}
    end
    arguments (Repeating)
        variables char
    end
    arguments
        options.Nt (1,1) double {mustBePositive} = Inf
        options.shouldOverwriteExisting double {mustBeMember(options.shouldOverwriteExisting,[0 1])} = 0 
        options.shouldAddDefaultVariables double {mustBeMember(options.shouldAddDefaultVariables,[0 1])} = 1 
        options.shouldUseClassicNetCDF double {mustBeMember(options.shouldUseClassicNetCDF,[0 1])} = 1 
    end


    [ncfile,matFilePath] = wvt.writeToFile(path,shouldOverwriteExisting=options.shouldOverwriteExisting,shouldAddDefaultVariables=0,shouldUseClassicNetCDF=options.shouldUseClassicNetCDF);

    % Now add a time dimension
    varAnnotation = wvt.variableAnnotationWithName('t');
    varAnnotation.attributes('units') = varAnnotation.units;
    varAnnotation.attributes('long_name') = varAnnotation.description;
    varAnnotation.attributes('standard_name') = 'time';
    varAnnotation.attributes('long_name') = 'time';
    varAnnotation.attributes('units') = 'seconds since 1970-01-01 00:00:00';
    varAnnotation.attributes('axis') = 'T';
    varAnnotation.attributes('calendar') = 'standard';
    ncfile.addDimension(varAnnotation.name,[],varAnnotation.attributes,options.Nt);

    if options.shouldAddDefaultVariables == 1
        if isempty(variables)
            variables = {'A0','Ap','Am'};
        else
            variables = union(variables,{'A0','Ap','Am'});
        end
    end

    for iVar=1:length(variables)
        if isKey(wvt.variableAnnotationNameMap,variables{iVar})
            varAnnotation = wvt.variableAnnotationWithName(variables{iVar});
        elseif isKey(wvt.propertyAnnotationNameMap,variables{iVar})
            varAnnotation = wvt.propertyAnnotationWithName(variables{iVar});
        else
            error('Unrecognized property or variable, %s',variables{iVar});
        end
        varAnnotation.attributes('units') = varAnnotation.units;
        varAnnotation.attributes('long_name') = varAnnotation.description;
        if varAnnotation.isComplex == 1
            ncfile.initComplexVariable(varAnnotation.name,horzcat(varAnnotation.dimensions,'t'),varAnnotation.attributes,'NC_DOUBLE');
        else
            ncfile.initVariable(varAnnotation.name,horzcat(varAnnotation.dimensions,'t'),varAnnotation.attributes,'NC_DOUBLE');
        end
    end

    outputIndex = 1;
    ncfile.concatenateVariableAlongDimension('t',wvt.t,'t',outputIndex);
    for iVar=1:length(variables)
        ncfile.concatenateVariableAlongDimension(variables{iVar},wvt.(variables{iVar}),'t',outputIndex);
    end

end