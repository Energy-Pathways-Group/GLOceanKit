function BuildDocumentationForWVTransformSubclass(classDocumentationFolder,parentName,parentFolder)
className = 'WVTransform';

targetFolder = sprintf('%s/%s',classDocumentationFolder,lower(className));
mc = meta.class.fromName(className);

metadataNameMap = ExtractMetadataFromClassPropertiesAndMethods(mc);

% Overwrite the automatically extracted metadata, with our version
dims = WVTransform.defaultDimensionAnnotations();
for iDim=1:length(dims)
    metadata = ExtractMetadataFromDetailedDescription(dims(iDim).detailedDescription);
    metadata.name = dims(iDim).name;
    metadata.className = className;
    metadata.shortDescription = dims(iDim).description;
    metadata.units = dims(iDim).units;
    metadata.isComplex = 0;
    metadata.functionType = FunctionType.transformDimension;
    metadataNameMap(metadata.name) = metadata;
end

% Overwrite the automatically extracted metadata, with our version
props = WVTransform.defaultPropertyAnnotations();
for iDim=1:length(props)
    metadata = ExtractMetadataFromDetailedDescription(props(iDim).detailedDescription);
    metadata.name = props(iDim).name;
    metadata.className = className;
    metadata.shortDescription = props(iDim).description;
    metadata.units = props(iDim).units;
    metadata.dimensions = props(iDim).dimensions;
    metadata.isComplex = props(iDim).isComplex;
    metadata.functionType = FunctionType.transformProperty;
    metadataNameMap(metadata.name) = metadata;
end

% Overwrite the automatically extracted metadata, with our version
ops = WVTransform.defaultOperations();
for iOp=1:length(ops)
    for iVar=1:length(ops(iOp).outputVariables)
        stateVar = ops(iOp).outputVariables(iVar);
        metadata = ExtractMetadataFromDetailedDescription(stateVar.detailedDescription);
        metadata.name = stateVar.name;
        metadata.className = className;
        metadata.shortDescription = stateVar.description;
        metadata.units = stateVar.units;
        metadata.dimensions = stateVar.dimensions;
        metadata.isComplex = stateVar.isComplex;
        metadata.functionType = FunctionType.stateVariable;
        metadataNameMap(metadata.name) = metadata;
    end
end

methods = WVTransform.defaultMethodAnnotations;
for iMethod=1:length(methods)
    metadata = ExtractMetadataFromDetailedDescription(methods(iMethod).detailedDescription);
    metadata.name = methods(iMethod).name;
    metadata.shortDescription = methods(iMethod).description;
    additionalMetadata = metadataNameMap(methods(iMethod).name);
    metadata.className = additionalMetadata.className;
    metadata.functionType = FunctionType.instanceMethod;
    metadataNameMap(metadata.name) = metadata;
end

[classDetailedDescription,classDefinedTopics] = ExtractTopicsFromClassAndSortMetadata(mc,metadataNameMap);

% Make a folder for all the class contents
if ~exist(targetFolder,'dir')
    mkdir(targetFolder);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Make the index/table of contents
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path = sprintf('%s/index.md',targetFolder);
MakeMarkdownFileForClass(path,className,classDetailedDescription,classDefinedTopics,metadataNameMap,parentName,parentFolder);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Make the individual pages
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iPageNumber = 0;
methodNames = keys(metadataNameMap);
for i=1:length(methodNames)
    iPageNumber = iPageNumber+1;
    metadata = metadataNameMap(methodNames{i});
    path = sprintf('%s/%s.md',targetFolder,lower(metadata.name));
    MakeMarkdownFileFromMethodMetadata(path,metadata,iPageNumber);
end


end