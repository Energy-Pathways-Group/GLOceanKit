cd(fileparts(which(mfilename)))

destinationFolder = '../../../../docs/';
sourceFolder = '../WebsiteDocumentation/';

copyfile(sourceFolder,destinationFolder);

parentName = 'WV transform & model';
parentFolder = 'classes-transform-and-model';
destinationFolder = '../../../../docs/classes-transform-and-model';

BuildDocumentationForClass(name='WVModel',folder=destinationFolder,parentName=parentName,parentFolder=parentFolder);
BuildDocumentationForWVTransformSubclass(destinationFolder,parentName,parentFolder);
BuildDocumentationForClass(name='WVTransformHydrostatic',folder=destinationFolder,parentName=parentName,parentFolder=parentFolder);
BuildDocumentationForClass(name='WVTransformConstantStratification',folder=destinationFolder,parentName=parentName,parentFolder=parentFolder);
BuildDocumentationForClass(name='WVTransformSingleMode',folder=destinationFolder,parentName=parentName,parentFolder=parentFolder);

parentName = 'Nonlinear flux operations';
parentFolder = 'classes-nonlinearfluxes';
destinationFolder = '../../../../docs/classes-nonlinearfluxes';

BuildDocumentationForClass(name='WVNonlinearFluxOperation',folder=destinationFolder,parentName=parentName,parentFolder=parentFolder);
BuildDocumentationForClass(name='Boussinesq',folder=destinationFolder,parentName=parentName,parentFolder=parentFolder);
BuildDocumentationForClass(name='BoussinesqSpatial',folder=destinationFolder,parentName=parentName,parentFolder=parentFolder);
BuildDocumentationForClass(name='QGPVE',folder=destinationFolder,parentName=parentName,parentFolder=parentFolder);

parentName = 'Operations & annotations';
parentFolder = 'classes-operations-and-annotations';
destinationFolder = '../../../../docs/classes-operations-and-annotations';

BuildDocumentationForClass(name='WVAnnotation',folder=destinationFolder,parentName=parentName,parentFolder=parentFolder);
BuildDocumentationForClass(name='WVDimensionAnnotation',folder=destinationFolder,parentName=parentName,parentFolder=parentFolder);
BuildDocumentationForClass(name='WVPropertyAnnotation',folder=destinationFolder,parentName=parentName,parentFolder=parentFolder);
BuildDocumentationForClass(name='WVVariableAnnotation',folder=destinationFolder,parentName=parentName,parentFolder=parentFolder);
BuildDocumentationForClass(name='WVOperation',folder=destinationFolder,parentName=parentName,parentFolder=parentFolder);

% ClassDocGenerator('WVFlowConstituent',classDocumentationFolder);



BuildDocumentationForClass(name='NetCDFFile',folder=destinationFolder,parentName=parentName,parentFolder=parentFolder);