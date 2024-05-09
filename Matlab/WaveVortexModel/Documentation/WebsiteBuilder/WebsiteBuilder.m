cd(fileparts(which(mfilename)))

buildFolder = '../../../../docs/';
sourceFolder = '../WebsiteDocumentation/';

copyfile(sourceFolder,buildFolder);

classFolderName = 'Class documentation';
websiteFolder = 'classes';
BuildDocumentationForWVTransformSubclass(buildFolder=buildFolder,websiteFolder=websiteFolder,parent=classFolderName,nav_order=1);
BuildDocumentationForClass(name='WVModel',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=classFolderName,nav_order=2);

parentName = 'Transforms';
websiteFolder = 'classes/transforms';

BuildDocumentationForClass(name='WVTransformBoussinesq',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName);
BuildDocumentationForClass(name='WVTransformHydrostatic',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName);
BuildDocumentationForClass(name='WVTransformConstantStratification',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName);
% BuildDocumentationForClass(name='WVTransformSingleMode',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName);

parentName = 'Nonlinear fluxes';
websiteFolder = 'classes/nonlinear-fluxes';

BuildDocumentationForClass(name='WVNonlinearFluxOperation',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName,nav_order=1);
BuildDocumentationForClass(name='WVNonlinearFlux',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName,nav_order=2);
BuildDocumentationForClass(name='WVNonlinearFluxForced',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName,nav_order=3);
BuildDocumentationForClass(name='WVNonlinearFluxQG',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName,nav_order=4);
BuildDocumentationForClass(name='WVNonlinearFluxQGForced',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName,nav_order=5);
BuildDocumentationForClass(name='WVNonlinearFluxWindForced',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName,nav_order=6);
BuildDocumentationForClass(name='WVNonlinearFluxSpatial',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName,nav_order=7);


parentName = 'Operations & annotations';
websiteFolder = 'classes/operations-and-annotations';

BuildDocumentationForClass(name='WVOperation',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName,nav_order=1);
BuildDocumentationForClass(name='WVAnnotation',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName,nav_order=2);
BuildDocumentationForClass(name='WVDimensionAnnotation',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName,nav_order=3);
BuildDocumentationForClass(name='WVPropertyAnnotation',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName,nav_order=4);
BuildDocumentationForClass(name='WVVariableAnnotation',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName,nav_order=5);

parentName = 'Flow components';
websiteFolder = 'classes/flow-components';

BuildDocumentationForClass(name='WVFlowComponent',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName,nav_order=1);
BuildDocumentationForClass(name='WVPrimaryFlowComponent',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName,nav_order=2);
BuildDocumentationForClass(name='WVGeostrophicComponent',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName,nav_order=3);
BuildDocumentationForClass(name='WVInternalGravityWaveComponent',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName,nav_order=4);
BuildDocumentationForClass(name='WVInertialOscillationComponent',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName,nav_order=5);
BuildDocumentationForClass(name='WVMeanDensityAnomalyComponent',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName,nav_order=6);

% ClassDocGenerator('WVFlowConstituent',classDocumentationFolder);



BuildDocumentationForClass(name='NetCDFFile',buildFolder=buildFolder,websiteFolder='classes',parent=classFolderName,nav_order=6);
BuildDocumentationForClass(name='WVGeometryDoublyPeriodic',buildFolder=buildFolder,websiteFolder='classes',parent=classFolderName,nav_order=7);