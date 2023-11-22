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
        options.shouldUseClassicNetCDF double {mustBeMember(options.shouldUseClassicNetCDF,[0 1])} = 1 
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
    ncfile = NetCDFFile(path,shouldUseClassicNetCDF=options.shouldUseClassicNetCDF);

    dims = {'x','y','z','k','l','j'};
    for iDim=1:length(dims)
        dimAnnotation = wvt.dimensionAnnotationWithName(dims{iDim});
        dimAnnotation.attributes('units') = dimAnnotation.units;
        dimAnnotation.attributes('long_name') = dimAnnotation.description;
        ncfile.addDimension(dimAnnotation.name,wvt.(dimAnnotation.name),dimAnnotation.attributes);
    end

    ncfile.addAttribute('source',sprintf('Created with the WaveVortexModel version %s',string(wvt.version)));
    ncfile.addAttribute('model_version',wvt.version);
    ncfile.addAttribute('date_created',string(datetime('now')));
    ncfile.addAttribute('history',string(strcat(string(datetime('now')),': file created.')));
    ncfile.addAttribute('references','Early, J., Lelong, M., & Sundermeyer, M. (2021). A generalized wave-vortex decomposition for rotating Boussinesq flows with arbitrary stratification. Journal of Fluid Mechanics, 912, A32. doi:10.1017/jfm.2020.995');
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
            varAnnotation = wvt.variableAnnotationWithName(variables{iVar});
        elseif isKey(wvt.propertyAnnotationNameMap,variables{iVar})
            varAnnotation = wvt.propertyAnnotationWithName(variables{iVar});
        else
            error('Unrecognized property or variable, %s',variables{iVar});
        end
        varAnnotation.attributes('units') = varAnnotation.units;
        varAnnotation.attributes('long_name') = varAnnotation.description;
        if varAnnotation.isComplex == 1
            ncfile.initComplexVariable(varAnnotation.name,varAnnotation.dimensions,varAnnotation.attributes,'NC_DOUBLE');
            ncfile.setVariable(varAnnotation.name,wvt.(varAnnotation.name));
        else
            ncfile.addVariable(varAnnotation.name,wvt.(varAnnotation.name),varAnnotation.dimensions,varAnnotation.attributes);
        end
    end

    ncfile.addAttribute('WVNonlinearFluxOperation',class(wvt.nonlinearFluxOperation));
    wvt.nonlinearFluxOperation.writeToFile(ncfile,wvt);
end