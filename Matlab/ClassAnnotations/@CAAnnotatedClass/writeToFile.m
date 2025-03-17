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

    if options.shouldAddRequiredProperties == 1
        properties = union(properties,self.requiredProperties);
    end

    propertyAnnotations = self.propertyAnnotationWithName(properties);

    ncfile = CAAnnotatedClass.writeToPath(self,path,shouldOverwriteExisting=options.shouldOverwriteExisting,propertyAnnotations=propertyAnnotations,attributes=options.attributes);
    ncfile.sync();
end