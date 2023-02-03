function [ncfile,matFilePath] = writeToFile(wvt,path,variables,options)
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
    arguments
        wvt WVTransform {mustBeNonempty}
        path char {mustBeNonempty}
    end
    arguments (Repeating)
        variables char
    end
    arguments
        options.shouldOverwriteExisting double {mustBeMember(options.shouldOverwriteExisting,[0 1])} = 0 
        options.shouldAddDefaultVariables double {mustBeMember(options.shouldAddDefaultVariables,[0 1])} = 1 
    end
    % Will not add 't' by default to allow for alternative definitions. Do
    % include 't' in the option input arguments if you want it written.
    [filepath,name,~] = fileparts(path);
    if isempty(filepath)
        matFilePath = sprintf('%s.mat',name);
    else
        matFilePath = sprintf('%s/%s.mat',filepath,name);
    end

    if options.shouldOverwriteExisting == 1
        if isfile(path)
            delete(path);
        end
        if isfile(matFilePath)
            delete(matFilePath);
        end
    else
        if isfile(path) || isfile(matFilePath)
            error('A file already exists with that name.')
        end
    end
    ncfile = NetCDFFile(path);

    dims = {'x','y','z','k','l','j'};
    for iDim=1:length(dims)
        transformDim = wvt.dimensionAnnotationWithName(dims{iDim});
        attributes = containers.Map();
        attributes('units') = transformDim.units;
        attributes('description') = transformDim.description;
        ncfile.addDimension(transformDim.name,wvt.(transformDim.name),attributes);
    end

    CreationDate = datestr(datetime('now'));
    ncfile.addAttribute('Model','Created with the WaveVortexModel, written by Jeffrey J. Early and collaborators.');
    ncfile.addAttribute('ModelVersion',wvt.version);
    ncfile.addAttribute('CreationDate',CreationDate);
    ncfile.addAttribute('WVTransform',class(wvt));

    attributesToWrite = {'latitude','t0','rho0','Lx','Ly','Lz'};
    variables = union(variables,attributesToWrite);

    if options.shouldAddDefaultVariables == 1
        if isempty(variables)
            variables = {'A0','Ap','Am','t'};
        else
            variables = union(variables,{'A0','Ap','Am','t'});
        end
    end

    for iVar=1:length(variables)
        if isKey(wvt.variableAnnotationNameMap,variables{iVar})
            transformVar = wvt.variableAnnotationWithName(variables{iVar});
        elseif isKey(wvt.propertyAnnotationNameMap,variables{iVar})
            transformVar = wvt.propertyAnnotationWithName(variables{iVar});
        else
            error('Unrecognized property or variable, %s',variables{iVar});
        end
        attributes = containers.Map();
        attributes('units') = transformVar.units;
        attributes('description') = transformVar.description;
        if transformVar.isComplex == 1
            ncfile.initComplexVariable(transformVar.name,transformVar.dimensions,attributes,'NC_DOUBLE');
            ncfile.setVariable(transformVar.name,wvt.(transformVar.name));
        else
            ncfile.addVariable(transformVar.name,wvt.(transformVar.name),transformVar.dimensions,attributes);
        end
    end
end