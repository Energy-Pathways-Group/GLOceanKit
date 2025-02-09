%%
atc_orig = AnnotationTestClass();
ncfile = atc_orig.writeToFile('test.nc',shouldOverwriteExisting=1);
atc_read = AnnotationTestClass.annotatedTestClassFromFile('test.nc');
assert(atc_orig.isequal(atc_read),"objects not equal");

%%
atc_orig = AnnotationTestClass();
atc_orig.myObjs(2) = AnnotationTestClassB();
ncfile = atc_orig.writeToFile('test.nc',shouldOverwriteExisting=1);
atc_read = AnnotationTestClass.annotatedTestClassFromFile('test.nc');
assert(atc_orig.isequal(atc_read),"objects not equal");

%%
atcB_orig = AnnotationTestClassB();
ncfile = atcB_orig.writeToFile('test.nc',shouldOverwriteExisting=1);
atcB_read = AnnotationTestClass.annotatedTestClassFromFile('test.nc');
assert(atcB_orig.isequal(atcB_read),"objects not equal");

%%
atc_orig = AnnotationTestClass();
newProp = CANumericProperty('newVar',{'x'},'', 'nothing');
atc_orig.addPropertyAnnotation(newProp);