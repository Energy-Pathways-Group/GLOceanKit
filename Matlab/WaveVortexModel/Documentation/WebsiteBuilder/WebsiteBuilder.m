cd(fileparts(which(mfilename)))

destinationFolder = '../../../../docs/';
sourceFolder = '../WebsiteDocumentation/';

copyfile(sourceFolder,destinationFolder);

return

transformAndModelClassDocumentationFolder = '../../../../docs/classes-transform-and-model';

BuildDocumentationForClass('WVModel',transformAndModelClassDocumentationFolder);
BuildDocumentationForWVTransformSubclass(transformAndModelClassDocumentationFolder);
BuildDocumentationForClass('WVTransformHydrostatic',transformAndModelClassDocumentationFolder);
BuildDocumentationForClass('WVTransformConstantStratification',transformAndModelClassDocumentationFolder);

nonlinearFluxClassDocumentationFolder = '../../../../docs/classes-nonlinearfluxes';

BuildDocumentationForClass('WVNonlinearFluxOperation',nonlinearFluxClassDocumentationFolder);
BuildDocumentationForClass('Boussinesq',nonlinearFluxClassDocumentationFolder);
BuildDocumentationForClass('BoussinesqSpatial',nonlinearFluxClassDocumentationFolder);

operationsClassDocumentationFolder = '../../../../docs/classes-operations-and-annotations';

BuildDocumentationForClass('WVAnnotation',operationsClassDocumentationFolder);
BuildDocumentationForClass('WVDimensionAnnotation',operationsClassDocumentationFolder);
BuildDocumentationForClass('WVPropertyAnnotation',operationsClassDocumentationFolder);
BuildDocumentationForClass('WVVariableAnnotation',operationsClassDocumentationFolder);
BuildDocumentationForClass('WVOperation',operationsClassDocumentationFolder);

% ClassDocGenerator('WVFlowConstituent',classDocumentationFolder);



BuildDocumentationForClass('NetCDFFile',operationsClassDocumentationFolder);