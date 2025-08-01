cd(fileparts(which(mfilename)))

buildFolder = '../../../../docs/';
sourceFolder = '../WebsiteDocumentation/';

copyfile(sourceFolder,buildFolder);

classFolderName = 'Class documentation';
websiteFolder = 'classes';
classDoc = WVTransformDocumentation('WVTransform',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=classFolderName,nav_order=1);
classDoc.writeToFile();

classDoc = ClassDocumentation('WVModel',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=classFolderName,nav_order=2);
classDoc.writeToFile();

%%
parentName = 'Transforms';
websiteFolder = 'classes/transforms';
classes = {'WVTransformBoussinesq','WVTransformHydrostatic','WVTransformConstantStratification'}; % 'WVTransformSingleMode'
excludedSuperclasses = {'handle','WVTransform'};
classDocumentation = WVTransformSubclassDocumentation.classDocumentationFromClassNames(classes,buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName,excludedSuperclasses=excludedSuperclasses);
arrayfun(@(a) a.writeToFile(),classDocumentation)

%%
parentName = 'Nonlinear fluxes';
websiteFolder = 'classes/nonlinear-fluxes';
classes = {'WVNonlinearFluxOperation','WVNonlinearFlux','WVNonlinearFluxForced','WVNonlinearFluxQG','WVNonlinearFluxQGForced','WVNonlinearFluxWindForced','WVNonlinearFluxSpatial'};

classDocumentation = ClassDocumentation.classDocumentationFromClassNames(classes,buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName);
arrayfun(@(a) a.writeToFile(),classDocumentation);

%%
parentName = 'Operations & annotations';
websiteFolder = 'classes/operations-and-annotations';
classes = {'WVOperation','WVAnnotation','WVDimensionAnnotation','WVPropertyAnnotation','WVVariableAnnotation'};

classDocumentation = ClassDocumentation.classDocumentationFromClassNames(classes,buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName);
arrayfun(@(a) a.writeToFile(),classDocumentation);


parentName = 'Flow components';
websiteFolder = 'classes/flow-components';
classes = {'WVFlowComponent','WVPrimaryFlowComponent','WVGeostrophicComponent','WVInternalGravityWaveComponent','WVInertialOscillationComponent','WVMeanDensityAnomalyComponent'};

classDocumentation = ClassDocumentation.classDocumentationFromClassNames(classes,buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName);
arrayfun(@(a) a.writeToFile(),classDocumentation);

classDoc = ClassDocumentation('NetCDFFile',buildFolder=buildFolder,websiteFolder='classes',parent=classFolderName,nav_order=6);
classDoc.writeToFile();

classDoc = ClassDocumentation('WVGeometryDoublyPeriodic',buildFolder=buildFolder,websiteFolder='classes',parent=classFolderName,nav_order=7);
classDoc.writeToFile();