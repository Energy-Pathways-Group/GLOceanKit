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

BuildDocumentationForClass(name='WVTransformHydrostatic',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName);
BuildDocumentationForClass(name='WVTransformConstantStratification',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName);
BuildDocumentationForClass(name='WVTransformSingleMode',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName);

parentName = 'Nonlinear fluxes';
websiteFolder = 'classes/nonlinear-fluxes';

BuildDocumentationForClass(name='WVNonlinearFluxOperation',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName,nav_order=1);
BuildDocumentationForClass(name='Boussinesq',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName,nav_order=2);
BuildDocumentationForClass(name='BoussinesqSpatial',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName,nav_order=3);
BuildDocumentationForClass(name='QGPVE',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName,nav_order=4);

parentName = 'Operations & annotations';
websiteFolder = 'classes/operations-and-annotations';

BuildDocumentationForClass(name='WVOperation',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName,nav_order=1);
BuildDocumentationForClass(name='WVAnnotation',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName,nav_order=2);
BuildDocumentationForClass(name='WVDimensionAnnotation',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName,nav_order=3);
BuildDocumentationForClass(name='WVPropertyAnnotation',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName,nav_order=4);
BuildDocumentationForClass(name='WVVariableAnnotation',buildFolder=buildFolder,websiteFolder=websiteFolder,parent=parentName,grandparent=classFolderName,nav_order=5);


% ClassDocGenerator('WVFlowConstituent',classDocumentationFolder);



BuildDocumentationForClass(name='NetCDFFile',buildFolder=buildFolder,websiteFolder='classes',parent=classFolderName,nav_order=6);