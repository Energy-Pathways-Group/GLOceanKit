cd(fileparts(which(mfilename)))

destinationFolder = '../../../../docs/';
sourceFolder = '../WebsiteDocumentation/';

copyfile(sourceFolder,destinationFolder);

classFolderName = 'Class documentation';
destinationFolder = '../../../../docs/classes';
BuildDocumentationForWVTransformSubclass(folder=destinationFolder,parent=classFolderName,nav_order=1);
BuildDocumentationForClass(name='WVModel',folder=destinationFolder,parent=classFolderName,nav_order=2);

parentName = 'Transforms';
destinationFolder = '../../../../docs/classes/transforms';

BuildDocumentationForClass(name='WVTransformHydrostatic',folder=destinationFolder,parent=parentName,grandparent=classFolderName);
BuildDocumentationForClass(name='WVTransformConstantStratification',folder=destinationFolder,parent=parentName,grandparent=classFolderName);
BuildDocumentationForClass(name='WVTransformSingleMode',folder=destinationFolder,parent=parentName,grandparent=classFolderName);

parentName = 'Nonlinear fluxes';
destinationFolder = '../../../../docs/classes/nonlinear-fluxes';

BuildDocumentationForClass(name='WVNonlinearFluxOperation',folder=destinationFolder,parent=parentName,grandparent=classFolderName,nav_order=1);
BuildDocumentationForClass(name='Boussinesq',folder=destinationFolder,parent=parentName,grandparent=classFolderName,nav_order=2);
BuildDocumentationForClass(name='BoussinesqSpatial',folder=destinationFolder,parent=parentName,grandparent=classFolderName,nav_order=3);
BuildDocumentationForClass(name='QGPVE',folder=destinationFolder,parent=parentName,grandparent=classFolderName,nav_order=4);

parentName = 'Operations & annotations';
destinationFolder = '../../../../docs/classes/operations-and-annotations';

BuildDocumentationForClass(name='WVOperation',folder=destinationFolder,parent=parentName,grandparent=classFolderName,nav_order=1);
BuildDocumentationForClass(name='WVAnnotation',folder=destinationFolder,parent=parentName,grandparent=classFolderName,nav_order=2);
BuildDocumentationForClass(name='WVDimensionAnnotation',folder=destinationFolder,parent=parentName,grandparent=classFolderName,nav_order=3);
BuildDocumentationForClass(name='WVPropertyAnnotation',folder=destinationFolder,parent=parentName,grandparent=classFolderName,nav_order=4);
BuildDocumentationForClass(name='WVVariableAnnotation',folder=destinationFolder,parent=parentName,grandparent=classFolderName,nav_order=5);


% ClassDocGenerator('WVFlowConstituent',classDocumentationFolder);



BuildDocumentationForClass(name='NetCDFFile',folder=destinationFolder,parent=classFolderName,nav_order=6);