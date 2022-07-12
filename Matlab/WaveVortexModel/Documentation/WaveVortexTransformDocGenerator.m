function WaveVortexTransformDocGenerator(classDocumentationFolder)
className = 'WaveVortexTransform';

targetFolder = sprintf('%s/%s',classDocumentationFolder,lower(className));
mc = meta.class.fromName(className);

detailedDescription = mc.DetailedDescription;

% extract topics and the detailed description (minus those topics)
subtopicExpression = '- topic:([ \t]*)(?<topicName>[^—\r\n]+)—([ \t]*)(?<subtopicName>[^\r\n]+)(?:$|\n)';
classDefinedSubTopics = regexpi(detailedDescription,subtopicExpression,'names');
detailedDescription = regexprep(detailedDescription,subtopicExpression,'','ignorecase');
topicExpression = '- topic:([ \t]*)(?<topicName>[^\r\n]+)(?:$|\n)';
classDefinedTopics = regexpi(detailedDescription,topicExpression,'names');
classDetailedDescription = regexprep(detailedDescription,topicExpression,'','ignorecase');

% classDefinedTopics is a structure array with field "topicName" and
% "subtopics".
% subtopics is itself a structure array with field "subtopicName" and
% "methodMetadata".
% methodMetadata is a structure array.
classDefinedTopicsIndex = containers.Map;
for iTopic=1:length(classDefinedTopics)
    classDefinedTopics(iTopic).topicName = strtrim(classDefinedTopics(iTopic).topicName);
    classDefinedTopicsIndex(classDefinedTopics(iTopic).topicName) = iTopic;
    classDefinedTopics(iTopic).subtopics = struct('subtopicName',[],'methodNames',{});
    classDefinedTopics(iTopic).subtopicsIndex = containers.Map;
end

for iSubtopic=1:length(classDefinedSubTopics)
    classDefinedSubTopics(iSubtopic).topicName = strtrim(classDefinedSubTopics(iSubtopic).topicName);
    classDefinedSubTopics(iSubtopic).subtopicName = strtrim(classDefinedSubTopics(iSubtopic).subtopicName);
    topicIndex = classDefinedTopicsIndex(classDefinedSubTopics(iSubtopic).topicName);
    classDefinedTopics(topicIndex).subtopics(end+1).subtopicName = classDefinedSubTopics(iSubtopic).subtopicName;
    classDefinedTopics(topicIndex).subtopicsIndex(classDefinedSubTopics(iSubtopic).subtopicName) = length(classDefinedTopics(topicIndex).subtopics);
end


methodsMap = containers.Map;
for i=1:length(mc.MethodList)
    metadata = ExtractMethodMetadata(mc.MethodList(i));
    if ~isempty(metadata)
        methodsMap(metadata.name) = metadata;
    end
end
for i=1:length(mc.PropertyList)
    metadata = ExtractMethodMetadata(mc.PropertyList(i));
    if ~isempty(metadata)
        methodsMap(metadata.name) = metadata;
    end
end

% Overwrite the automatically extracted metadata, with our version
dims = WaveVortexTransform.defaultTransformDimensions();
for iDim=1:length(dims)
    metadata = ExtractMetadataFromDetailedDescription(dims(iDim).detailedDescription);
    metadata.name = dims(iDim).name;
    metadata.className = className;
    metadata.shortDescription = dims(iDim).description;
    metadata.units = dims(iDim).units;
    methodsMap(metadata.name) = metadata;
end

% Overwrite the automatically extracted metadata, with our version
props = WaveVortexTransform.defaultTransformProperties();
for iDim=1:length(props)
    metadata = ExtractMetadataFromDetailedDescription(props(iDim).detailedDescription);
    metadata.name = props(iDim).name;
    metadata.className = className;
    metadata.shortDescription = props(iDim).description;
    metadata.units = props(iDim).units;
    metadata.dimensions = props(iDim).dimensions;
    metadata.isComplex = props(iDim).isComplex;
    methodsMap(metadata.name) = metadata;
end

% Overwrite the automatically extracted metadata, with our version
ops = WaveVortexTransform.defaultTransformOperations();
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
        methodsMap(metadata.name) = metadata;
    end
end

methodNames = keys(methodsMap);
for i=1:length(methodNames)
    metadata = methodsMap(methodNames{i});

    % Force everything to have a subtopic
    if ~isfield(metadata,'topic') || isempty(metadata.topic)
        metadata.topic = 'Other';
    end

    if ~isfield(metadata,'subtopic') || isempty(metadata.subtopic)
        metadata.subtopic = 'Other';
    end

    % If this topic isn't created, create it!
    if ~isKey(classDefinedTopicsIndex,metadata.topic)
        classDefinedTopics(end+1).topicName = metadata.topic;
        topicIndex = length(classDefinedTopics);
        classDefinedTopicsIndex(metadata.topic) = topicIndex;
        classDefinedTopics(topicIndex).subtopics = struct('subtopicName',[],'methodNames',{});
        classDefinedTopics(topicIndex).subtopicsIndex = containers.Map;
    end
    topicIndex = classDefinedTopicsIndex(metadata.topic);

    % If the subtopic isn't created, create it!
    if ~isKey(classDefinedTopics(topicIndex).subtopicsIndex,metadata.subtopic)
        classDefinedTopics(topicIndex).subtopics(end+1).subtopicName = metadata.subtopic;
        classDefinedTopics(topicIndex).subtopicsIndex(metadata.subtopic) = length(classDefinedTopics(topicIndex).subtopics);
    end
    subtopicIndex = classDefinedTopics(topicIndex).subtopicsIndex(metadata.subtopic);
    
    % And finally, add this method name to the appropriate subtopic index
    classDefinedTopics(topicIndex).subtopics(subtopicIndex).methodNames{end+1} = metadata.name;
end

% Make a folder for all the class contents
if ~exist(targetFolder,'dir')
    mkdir(targetFolder);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Make the index/table of contents
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileID = fopen(sprintf('%s/index.md',targetFolder),'w');
fprintf(fileID,'---\nlayout: default\ntitle: %s\nparent: Classes\nhas_children: false\nhas_toc: false\n---\n\n',className);
fprintf(fileID,'#  %s\n',className);
fprintf(fileID,'\n%s\n\n',mc.Description);
if ~isempty(mc.DetailedDescription)
    fprintf(fileID,'## Discussion\n%s\n',classDetailedDescription);
end
fprintf(fileID,'\n\n## Topics\n');

for topicIndex = 1:length(classDefinedTopics)
    fprintf(fileID,'+ %s\n',classDefinedTopics(topicIndex).topicName);

    % If the 'Other' subtopic isn't there, it means no methods actually
    % used this topic.
    if ~isKey(classDefinedTopics(topicIndex).subtopicsIndex,'Other')
        continue;
    end

    otherSubtopicIndex = classDefinedTopics(topicIndex).subtopicsIndex('Other');
    for subtopicIndex = 1:length(classDefinedTopics(topicIndex).subtopics)
        subtopic = classDefinedTopics(topicIndex).subtopics(subtopicIndex);
        fprintf(fileID,'  + %s\n',subtopic.subtopicName);
        for methodIndex = 1:length(subtopic.methodNames)
            fprintf(fileID,'    + [`%s`](/classes/%s/%s.html) ',subtopic.methodNames{methodIndex},lower(className),lower(subtopic.methodNames{methodIndex}));
            fprintf(fileID,'%s\n',methodsMap(methodNames{methodIndex}).shortDescription);
        end
    end
end

fprintf(fileID,'\n\n---');

fclose(fileID);



end