function atc = annotatedClassFromFile(path)
ncfile = NetCDFFile(path);
if isKey(ncfile.attributes,'AnnotatedClass')
    className = ncfile.attributes('AnnotatedClass');
    requiredProperties = feval(strcat(className,'.classRequiredPropertyNames'));
    for iVar = 1:length(requiredProperties)
        name = requiredProperties{iVar};
        var.(name) = ncfile.readVariables(name);
    end
    varCell = namedargs2cell(var);
    atc = feval(className,varCell{:});
else
    error('Unable to find the attribute AnnotatedClass');
end
end