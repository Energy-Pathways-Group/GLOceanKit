cd(fileparts(which(mfilename)))

destinationFolder = '../../../../docs/';
sourceFolder = '../WebsiteDocumentation/';

copyfile(sourceFolder,destinationFolder);

folderName = 'WV transform & model';
destinationFolder = '../../../../docs/classes-transform-and-model';

BuildDocumentationForClass(name='WVModel',folder=destinationFolder,parent=folderName);
BuildDocumentationForWVTransformSubclass(destinationFolder,folderName);
BuildDocumentationForClass(name='WVTransformHydrostatic',folder=destinationFolder,parent=folderName);
BuildDocumentationForClass(name='WVTransformConstantStratification',folder=destinationFolder,parent=folderName);
BuildDocumentationForClass(name='WVTransformSingleMode',folder=destinationFolder,parent=folderName);

folderName = 'Nonlinear flux operations';
destinationFolder = '../../../../docs/classes-nonlinearfluxes';

BuildDocumentationForClass(name='WVNonlinearFluxOperation',folder=destinationFolder,parent=folderName);
BuildDocumentationForClass(name='Boussinesq',folder=destinationFolder,parent=folderName);
BuildDocumentationForClass(name='BoussinesqSpatial',folder=destinationFolder,parent=folderName);
BuildDocumentationForClass(name='QGPVE',folder=destinationFolder,parent=folderName);

folderName = 'Operations & annotations';
destinationFolder = '../../../../docs/classes-operations-and-annotations';

BuildDocumentationForClass(name='WVAnnotation',folder=destinationFolder,parent=folderName);
BuildDocumentationForClass(name='WVDimensionAnnotation',folder=destinationFolder,parent=folderName);
BuildDocumentationForClass(name='WVPropertyAnnotation',folder=destinationFolder,parent=folderName);
BuildDocumentationForClass(name='WVVariableAnnotation',folder=destinationFolder,parent=folderName);
BuildDocumentationForClass(name='WVOperation',folder=destinationFolder,parent=folderName);

% ClassDocGenerator('WVFlowConstituent',classDocumentationFolder);



BuildDocumentationForClass(name='NetCDFFile',folder=destinationFolder,parent=folderName);