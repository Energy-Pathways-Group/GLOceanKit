function ncfile = writeToFile(self,path,properties,options)
    % Write this instance to NetCDF file.
    %
    % Writes the WVTransform instance to file, with enough information to
    % re-initialize. Pass additional variables to the variable list that
    % should also be written to file.
    %
    % Subclasses should add any necessary properties or variables to the
    % variable list before calling this superclass method.
    %
    % - Topic: Write to file
    % - Declaration: ncfile = writeToFile(path,properties,options)
    % - Parameter path: path to write file
    % - Parameter properties: strings of property names.
    % - Parameter shouldOverwriteExisting: (optional) boolean indicating whether or not to overwrite an existing file at the path. Default 0. 
    % - Parameter shouldAddDefaultVariables: (optional) boolean indicating whether or not add default variables `A0`,`Ap`,`Am`,`t`. Default 1.
    % - Parameter attributes: (optional) dictionary containing additional attributes to add thet NetCDF file
    arguments (Input)
        self CAAnnotatedClass {mustBeNonempty}
        path char {mustBeNonempty}
    end
    arguments (Input,Repeating)
        properties char
    end
    arguments (Input)
        options.shouldOverwriteExisting logical = false
        options.shouldAddRequiredProperties logical = true
        options.attributes = configureDictionary("string","string")
    end
    arguments (Output)
        ncfile NetCDFFile
    end
    options.attributes('WVTransform') = class(self);
    options.attributes('source') = sprintf('Created with the WaveVortexModel version %s',string(self.version));
    options.attributes('model_version') = self.version;
    options.attributes('date_created') = string(datetime('now'));
    options.attributes('history') = string(strcat(string(datetime('now')),': file created.'));
    options.attributes('references') = 'Early, J., Lelong, M., & Sundermeyer, M. (2021). A generalized wave-vortex decomposition for rotating Boussinesq flows with arbitrary stratification. Journal of Fluid Mechanics, 912, A32. doi:10.1017/jfm.2020.995';

    optionCell = namedargs2cell(options);
    ncfile = writeToFile@CAAnnotatedClass(self,path,properties{:},optionCell{:});
end