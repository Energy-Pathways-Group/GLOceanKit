function [detailedDescription,rootTopic] = ExtractTopicsFromClassAndSortMetadata(mc,metadataNameMap)
% Extracts topics and subtopic from a detailedDescription and creates a structure useful creating an ordered topic index
%
% - Parameter mc: the detailed description 
% - Parameter metadataNameMap: map with keys of method names and values of normalized metadata structs. 
% - Returns detailedDescription: the detailed description without the topic/subtopic metadata
% - Returns classDefinedTopics: a structure array with field "topicName" (char), "subtopics" (struct('subtopicName',[],'methodNames',{})), and "subtopicsIndex" (containers.Map with subtopic name mapping to an order index).
arguments
    mc meta.class
    metadataNameMap containers.Map
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Step 1) extract a list of topics and subtopics from the class's
% detailedDescription, and then sort those into useful structures.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

detailedDescription = mc.DetailedDescription;

% extract topics and the detailed description (minus those topics)
subsubtopicExpression = '- topic:([ \t]*)(?<topicName>[^—\r\n]+)—([ \t]*)(?<subtopicName>[^\r\n]+)—([ \t]*)(?<subsubtopicName>[^\r\n]+)(?:$|\n)';
classDefinedSubsubTopics = regexpi(detailedDescription,subsubtopicExpression,'names');
detailedDescription = regexprep(detailedDescription,subsubtopicExpression,'','ignorecase');
subtopicExpression = '- topic:([ \t]*)(?<topicName>[^—\r\n]+)—([ \t]*)(?<subtopicName>[^\r\n]+)(?:$|\n)';
classDefinedSubTopics = regexpi(detailedDescription,subtopicExpression,'names');
detailedDescription = regexprep(detailedDescription,subtopicExpression,'','ignorecase');
topicExpression = '- topic:([ \t]*)(?<topicName>[^\r\n]+)(?:$|\n)';
classDefinedTopics = regexpi(detailedDescription,topicExpression,'names');
detailedDescription = regexprep(detailedDescription,topicExpression,'','ignorecase');

rootTopic = Topic('Root');
for iTopic=1:length(classDefinedTopics)
    rootTopic.addSubtopic(Topic(strtrim(classDefinedTopics(iTopic).topicName)));
end

for iSubtopic=1:length(classDefinedSubTopics)
    topic = rootTopic.subtopicWithName(strtrim(classDefinedSubTopics(iSubtopic).topicName));
    topic.addSubtopic(Topic(strtrim(classDefinedSubTopics(iSubtopic).subtopicName)));
end

for iSubtopic=1:length(classDefinedSubsubTopics)
    topic = rootTopic.subtopicWithName(strtrim(classDefinedSubsubTopics(iSubtopic).topicName));
    subtopic = topic.subtopicWithName(strtrim(classDefinedSubsubTopics(iSubtopic).subtopicName));
    subtopic.addSubtopic(Topic(strtrim(classDefinedSubsubTopics(iSubtopic).subsubtopicName)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Step 2) take the method/property/dimension/etc metadata sort those into
% the appropriate topics.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

methodNames = keys(metadataNameMap);
otherTopic = Topic('Other');
for i=1:length(methodNames)
    metadata = metadataNameMap(methodNames{i});

    % if the method has no topic, assign it to the 'other' topic we created
    if ~isfield(metadata,'topic') || isempty(metadata.topic)
        otherTopic.addMethod(metadata);
        continue;
    end

    % if the topic does not exist, create it
    if isempty(rootTopic.subtopicWithName(metadata.topic))
        rootTopic.addSubtopic(Topic(metadata.topic));
    end
    topic = rootTopic.subtopicWithName(metadata.topic);

    % if there is no subtopic, assign to the topic, and then exit
    if ~isfield(metadata,'subtopic') || isempty(metadata.subtopic)
        topic.addMethod(metadata);
        continue;
    end

    % if the subtopic does not exist, create it
    if isempty(topic.subtopicWithName(metadata.subtopic))
        topic.addSubtopic(Topic(metadata.subtopic));
    end
    subtopic = topic.subtopicWithName(metadata.subtopic);

    % if there is no subsubtopic, assign to the subtopic, and then exit
    if ~isfield(metadata,'subsubtopic') || isempty(metadata.subsubtopic)
        subtopic.addMethod(metadata);
        continue;
    end
    
    % if the subsubtopic does not exist, create it
    if isempty(subtopic.subtopicWithName(metadata.subsubtopic))
        subtopic.addSubtopic(Topic(metadata.subsubtopic));
    end
    subsubtopic = subtopic.subtopicWithName(metadata.subsubtopic);
    subsubtopic.addMethod(metadata);
end

% Make the 'Other' topic go at the end
rootTopic.addSubtopic(otherTopic);