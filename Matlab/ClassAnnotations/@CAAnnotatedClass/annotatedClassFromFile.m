function atc = annotatedClassFromFile(path)
ncfile = NetCDFFile(path);
if isKey(ncfile.attributes,'AnnotatedClass')
    className = ncfile.attributes('AnnotatedClass');
    if ncfile.hasGroupWithName(className)
        group = ncfile.groupWithName(className);
    else
        group = ncfile;
    end
    atc = CAAnnotatedClass.annotatedClassFromGroup(group);
else
    error('Unable to find the attribute AnnotatedClass');
end
end