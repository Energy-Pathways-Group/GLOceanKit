function BuildDocumentationForClass(options)
% Generates documentation for a class from the class metadata
%
% Specifically, this script takes the class metadata and creates a markdown
% page for that class using the class's detailed description as the primary
% text, followed by an index, sorted by topic, to all the properties and
% methods of the class.
%
% The metadata for the property and methods is extracted from the property
% and method metadata (from Matlab's metaclass) as well. Each property and
% method gets its own page.
arguments
    options.name
    options.folder
    options.parentName
    options.parentFolder
end

mc = meta.class.fromName(options.name);

metadataNameMap = ExtractMetadataFromClassPropertiesAndMethods(mc);

[classDetailedDescription,classDefinedTopics] = ExtractTopicsFromClassAndSortMetadata(mc,metadataNameMap);


% Make a folder for all the class contents
targetFolder = sprintf('%s/%s',options.folder,lower(options.name));
if ~exist(targetFolder,'dir')
    mkdir(targetFolder);
end
path = sprintf('%s/index.md',targetFolder);


MakeMarkdownFileForClass(path,options.name,classDetailedDescription,classDefinedTopics,metadataNameMap,options.parentName,options.parentFolder);


iPageNumber = 0;
methodNames = keys(metadataNameMap);
for i=1:length(methodNames)
    iPageNumber = iPageNumber+1;
    metadata = metadataNameMap(methodNames{i});
    path = sprintf('%s/%s.md',targetFolder,lower(metadata.name));
    MakeMarkdownFileFromMethodMetadata(path,metadata,iPageNumber);
end

end