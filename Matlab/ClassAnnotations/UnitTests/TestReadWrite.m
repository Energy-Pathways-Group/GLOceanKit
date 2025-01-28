%%
a = AnnotationTestClass();
ncfile = a.writeToFile('test.nc',shouldOverwriteExisting=1);
b = AnnotationTestClass.annotatedTestClassFromFile('test.nc');

%%
a = AnnotationTestClass();
newProp = CANumericProperty('newVar',{'x'},'', 'nothing');
a.addPropertyAnnotation(newProp);