cd(fileparts(which(mfilename)))
classDocumentationFolder = '../../../../docs/classes';

BuildDocumentationForClass('WVAnnotation',classDocumentationFolder);
BuildDocumentationForClass('WVDimensionAnnotation',classDocumentationFolder);
BuildDocumentationForClass('WVPropertyAnnotation',classDocumentationFolder);
BuildDocumentationForClass('WVVariableAnnotation',classDocumentationFolder);
BuildDocumentationForClass('WVOperation',classDocumentationFolder);
BuildDocumentationForClass('WVNonlinearFluxOperation',classDocumentationFolder);

% ClassDocGenerator('WVFlowConstituent',classDocumentationFolder);

BuildDocumentationForClass('WVModel',classDocumentationFolder);

BuildDocumentationForWVTransformSubclass(classDocumentationFolder);
BuildDocumentationForClass('WVTransformHydrostatic',classDocumentationFolder);
BuildDocumentationForClass('WVTransformConstantStratification',classDocumentationFolder);

BuildDocumentationForClass('NetCDFFile',classDocumentationFolder);