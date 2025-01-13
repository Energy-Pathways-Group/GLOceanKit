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
    arguments (Input)
        wvt WVTransform {mustBeNonempty}
        path char {mustBeNonempty}
    end
    arguments (Input,Repeating)
        variables char
    end
    arguments (Input)
        options.shouldOverwriteExisting double {mustBeMember(options.shouldOverwriteExisting,[0 1])} = 0 
        options.shouldAddDefaultVariables double {mustBeMember(options.shouldAddDefaultVariables,[0 1])} = 1 
        options.shouldUseClassicNetCDF double {mustBeMember(options.shouldUseClassicNetCDF,[0 1])} = 1 
        options.dimensions = {}
    end
    arguments (Output)
        ncfile NetCDFFile
        matFilePath char
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

    if ~iscell(options.dimensions)
        if isempty(options.dimensions)
            dimensions = {};
        else
            dimensions = {options.dimensions};
        end
    else
        dimensions = options.dimensions;
    end

    dims = union(dimensions,{'x','y','z','kl','j'});
    for iDim=1:length(dims)
        dimAnnotation = wvt.dimensionAnnotationWithName(dims{iDim});
        dimAnnotation.attributes('units') = dimAnnotation.units;
        dimAnnotation.attributes('long_name') = dimAnnotation.description;
        ncfile.addDimension(dimAnnotation.name,wvt.(dimAnnotation.name),attributes=dimAnnotation.attributes);
    end

    ncfile.addAttribute('source',sprintf('Created with the WaveVortexModel version %s',string(wvt.version)));
    ncfile.addAttribute('model_version',wvt.version);
    ncfile.addAttribute('date_created',string(datetime('now')));
    ncfile.addAttribute('history',string(strcat(string(datetime('now')),': file created.')));
    ncfile.addAttribute('references','Early, J., Lelong, M., & Sundermeyer, M. (2021). A generalized wave-vortex decomposition for rotating Boussinesq flows with arbitrary stratification. Journal of Fluid Mechanics, 912, A32. doi:10.1017/jfm.2020.995');
    ncfile.addAttribute('WVTransform',class(wvt));

    % Pull out the user requested attributes that should be written to the
    % global attributes of the ncfile
    newAttributes = {};
    for iVar=1:length(variables)
        if isKey(wvt.propertyAnnotationNameMap,variables{iVar})
            newAttributes{end+1} = variables{iVar};
        end
    end
    variables = setdiff(variables,newAttributes);

    % Write these attributes to the root group
    attributesToWrite = union(newAttributes,{'latitude','t0','rho0','Lx','Ly','Lz','k','l','shouldAntialias'});
    for iVar=1:length(attributesToWrite)
        varAnnotation = wvt.propertyAnnotationWithName(attributesToWrite{iVar});
        varAnnotation.attributes('units') = varAnnotation.units;
        varAnnotation.attributes('long_name') = varAnnotation.description;
        ncfile.addVariable(varAnnotation.name,varAnnotation.dimensions,wvt.(varAnnotation.name),isComplex=varAnnotation.isComplex,attributes=varAnnotation.attributes);
    end

    % Any variables should be written with a possible time dimension to the
    % "wv" group.
    if options.shouldAddDefaultVariables == 1
        if isempty(variables)
            variables = {'A0','Ap','Am','t'};
        else
            variables = union(variables,{'A0','Ap','Am','t'});
        end
    end

    group = ncfile.addGroup("wv");
    for iVar=1:length(variables)
        varAnnotation = wvt.variableAnnotationWithName(variables{iVar});
        varAnnotation.attributes('units') = varAnnotation.units;
        varAnnotation.attributes('long_name') = varAnnotation.description;
        group.addVariable(varAnnotation.name,varAnnotation.dimensions,wvt.(varAnnotation.name),isComplex=varAnnotation.isComplex,attributes=varAnnotation.attributes);
    end

    ncfile.addAttribute("TotalForcingGroups",length(self.forcing))
    for iForce=1:length(self.forcing)
        forceGroup = ncfile.addGroup("forcing-"+iForce);
        forceGroup.addAttribute('WVForcing',class(self.forcing{iForce}));
        self.forcing{iForce}.writeToFile(forceGroup,wvt);
    end
end