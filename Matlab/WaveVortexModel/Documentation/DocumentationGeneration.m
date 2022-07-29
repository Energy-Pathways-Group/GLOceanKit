cd(fileparts(which(mfilename)))
classDocumentationFolder = '../../../docs/classes';

ClassDocGenerator('WVAnnotation',classDocumentationFolder);
ClassDocGenerator('WVDimensionAnnotation',classDocumentationFolder);
ClassDocGenerator('WVPropertyAnnotation',classDocumentationFolder);
ClassDocGenerator('WVVariableAnnotation',classDocumentationFolder);

ClassDocGenerator('WVOperation',classDocumentationFolder);
ClassDocGenerator('WVNonlinearFluxOperation',classDocumentationFolder);

% ClassDocGenerator('WVFlowConstituents',classDocumentationFolder);

ClassDocGenerator('WVModel',classDocumentationFolder);

WaveVortexTransformDocGenerator(classDocumentationFolder);