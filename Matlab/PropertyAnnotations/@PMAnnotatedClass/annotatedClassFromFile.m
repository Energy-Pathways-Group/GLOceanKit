function atc = annotatedClassFromFile(path)
ncfile = NetCDFFile(path);
if isKey(ncfile.attributes,'PMAnnotatedClass')
    className = ncfile.attributes('PMAnnotatedClass');
    requiredVariables = union(feval(strcat(className,'.classRequiredDimensions')),feval(strcat(className,'.classRequiredVariables')));
    for iVar = 1:length(requiredVariables)
        name = requiredVariables{iVar};
        var.(name) = ncfile.readVariables(name);
    end
    varCell = namedargs2cell(var);
    atc = feval(className,varCell{:});
else
    error('Unable to find the attribute PMAnnotatedClass');
end
end