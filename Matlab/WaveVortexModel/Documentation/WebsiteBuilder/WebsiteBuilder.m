cd(fileparts(which(mfilename)))

transformAndModelClassDocumentationFolder = '../../../../docs/classestransformandmodel';

BuildDocumentationForClass('WVModel',classDocumentationFolder);
BuildDocumentationForWVTransformSubclass(classDocumentationFolder);
BuildDocumentationForClass('WVTransformHydrostatic',classDocumentationFolder);
BuildDocumentationForClass('WVTransformConstantStratification',classDocumentationFolder);

nonlinearFluxClassDocumentationFolder = '../../../../docs/classesnonlinearflux';

BuildDocumentationForClass('WVNonlinearFluxOperation',classDocumentationFolder);
BuildDocumentationForClass('Boussinesq',classDocumentationFolder);
BuildDocumentationForClass('BoussinesqSpatial',classDocumentationFolder);

classDocumentationFolder = '../../../../docs/classes';

BuildDocumentationForClass('WVAnnotation',classDocumentationFolder);
BuildDocumentationForClass('WVDimensionAnnotation',classDocumentationFolder);
BuildDocumentationForClass('WVPropertyAnnotation',classDocumentationFolder);
BuildDocumentationForClass('WVVariableAnnotation',classDocumentationFolder);
BuildDocumentationForClass('WVOperation',classDocumentationFolder);
BuildDocumentationForClass('WVNonlinearFluxOperation',classDocumentationFolder);

% ClassDocGenerator('WVFlowConstituent',classDocumentationFolder);



BuildDocumentationForClass('NetCDFFile',classDocumentationFolder);