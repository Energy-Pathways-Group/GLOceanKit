a = AnnotationTestClass();
ncfile = a.writeToFile('test.nc',shouldOverwriteExisting=1);
b = AnnotationTestClass.annotatedTestClassFromFile('test.nc');