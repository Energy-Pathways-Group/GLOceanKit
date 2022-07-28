cd(fileparts(which(mfilename)))

className = 'WVModel';
classDocumentationFolder = '../../../docs/classes';

ClassDocGenerator(className,classDocumentationFolder);

WaveVortexTransformDocGenerator(classDocumentationFolder);