cd(fileparts(which(mfilename)))
classDocumentationFolder = '../../../docs/classes';

ClassDocGenerator('WVAnnotation',classDocumentationFolder);
ClassDocGenerator('WVDimensionAnnotation',classDocumentationFolder);
ClassDocGenerator('WVPropertyAnnotation',classDocumentationFolder);

% ClassDocGenerator('WVFlowConstituents',classDocumentationFolder);

ClassDocGenerator('WVModel',classDocumentationFolder);

WaveVortexTransformDocGenerator(classDocumentationFolder);