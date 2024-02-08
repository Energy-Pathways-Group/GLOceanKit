function concatenateVariablesAlongTimeDimension(wvt,path)
    % Concatenate variables along the time dimension
    %
    % Writes to an existing file initialized with
    % createNetCDFFileForTimeStepOutput.
    %
    % - Topic: Write to file
    % - Declaration: concatenateVariablesAlongTimeDimension(path)
    % - Parameter path: path to write file
    arguments
        wvt WVTransform {mustBeNonempty}
        path char {mustBeNonempty}
    end
    ncfile = NetCDFFile(path);
    for iDim = 1:length(ncfile.dimensions)
        if strcmp(ncfile.dimensions(iDim).name,'t')
            outputIndex = ncfile.dimensions(iDim).nPoints + 1;
        end
    end

    timeDependentVariables = {};
    for iVar = 1:length(ncfile.complexVariables)
        for iDim = 1:length(ncfile.complexVariables(iVar).realVar.dimensions)
            if strcmp(ncfile.complexVariables(iVar).realVar.dimensions(iDim).name,'t')
                timeDependentVariables = union(timeDependentVariables,ncfile.complexVariables(iVar).name);
            end
        end
    end

    for iVar = 1:length(ncfile.variables)
        if isKey(ncfile.variables(iVar).attributes,{ncfile.GLNetCDFSchemaIsComplexKey})
            continue;
        end
        if strcmp(ncfile.variables(iVar).name,'t')
            continue;
        end

        for iDim = 1:length(ncfile.variables(iVar).dimensions)
            if strcmp(ncfile.variables(iVar).dimensions(iDim).name,'t')
                timeDependentVariables = union(timeDependentVariables,ncfile.variables(iVar).name);
            end
        end
    end

    ncfile.concatenateVariableAlongDimension('t',wvt.t,'t',outputIndex);
    for iVar=1:length(timeDependentVariables)
        ncfile.concatenateVariableAlongDimension(timeDependentVariables{iVar},wvt.(timeDependentVariables{iVar}),'t',outputIndex);
    end
end