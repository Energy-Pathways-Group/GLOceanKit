function ncfile = writeToFile(wvt,path,variables,options)
    % Output the `WaveVortexTransform` to file.
    %
    % - Topic: Write to file
    % - Declaration: ncfile = writeToFile(netcdfFile,variables,options)
    % - Parameter path: path to write file
    % - Parameter variables: strings of variable names.
    % - Parameter shouldOverwriteExisting: (optional) boolean indicating whether or not to overwrite an existing file at the path. Default 0. 
    % - Parameter shouldAddDefaultVariables: (optional) boolean indicating whether or not add default variables `A0`,`Ap`,`Am`,`t`. Default 1.
    arguments
        wvt WaveVortexTransform {mustBeNonempty}
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
    ncfile.addAttribute('WaveVortexTransform',class(wvt));

    attributesToWrite = {'latitude','t0','rho0','Lx','Ly','Lz'};

    if isa(wvt,'WaveVortexTransformConstantStratification')
        attributesToWrite = union({'N0'},attributesToWrite);
    elseif isa(wvt,'WaveVortexTransformSingleMode')
        attributesToWrite = union({'h'},attributesToWrite);
    elseif isa(wvt,'WaveVortexTransformHydrostatic')
        %                 attributesToWrite = union({'rhobar','N2','dLnN2','PFinv','QGinv','PF','QG','h','P','Q'},attributesToWrite);
        attributesToWrite = union({'rhobar','N2','dLnN2','PFinv','QGinv','PF','QG'},attributesToWrite);

        rhoFunction = wvt.rhoFunction;
        N2Function = wvt.N2Function;
        dLnN2Function = wvt.dLnN2Function;
        save(matFilePath,'rhoFunction','N2Function','dLnN2Function','CreationDate');
        fprintf('In addition to the NetCDF file, a .mat sidecar file was created at the same path.\n');
    else
        error('Not implemented');
    end

    for iVar=1:length(attributesToWrite)
        transformVar = wvt.propertyAnnotationWithName(attributesToWrite{iVar});
        attributes = containers.Map();
        attributes('units') = transformVar.units;
        attributes('description') = transformVar.description;
        ncfile.addVariable(transformVar.name,wvt.(transformVar.name),transformVar.dimensions,attributes);
    end

    if options.shouldAddDefaultVariables == 1
        if isempty(variables)
            variables = {'A0','Ap','Am','t'};
        else
            variables = union(variables,{'A0','Ap','Am','t'});
        end
    end

    for iVar=1:length(variables)
        transformVar = wvt.variableAnnotationWithName(variables{iVar});
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