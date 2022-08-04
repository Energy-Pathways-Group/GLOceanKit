function [detailedDescription,classDefinedTopics] = ExtractTopicsFromClassAndSortMetadata(mc,metadataNameMap)
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
subtopicExpression = '- topic:([ \t]*)(?<topicName>[^—\r\n]+)—([ \t]*)(?<subtopicName>[^\r\n]+)(?:$|\n)';
classDefinedSubTopics = regexpi(detailedDescription,subtopicExpression,'names');
detailedDescription = regexprep(detailedDescription,subtopicExpression,'','ignorecase');
topicExpression = '- topic:([ \t]*)(?<topicName>[^\r\n]+)(?:$|\n)';
classDefinedTopics = regexpi(detailedDescription,topicExpression,'names');
detailedDescription = regexprep(detailedDescription,topicExpression,'','ignorecase');

% classDefinedTopicsIndex: containers.Map with key the topic name and value the topic order name
classDefinedTopicsIndex = containers.Map;
for iTopic=1:length(classDefinedTopics)
    classDefinedTopics(iTopic).topicName = strtrim(classDefinedTopics(iTopic).topicName);
    classDefinedTopicsIndex(lower(classDefinedTopics(iTopic).topicName)) = iTopic;
    classDefinedTopics(iTopic).subtopics = struct('subtopicName',[],'methodNames',{});
    classDefinedTopics(iTopic).subtopicsIndex = containers.Map;
end

for iSubtopic=1:length(classDefinedSubTopics)
    classDefinedSubTopics(iSubtopic).topicName = strtrim(classDefinedSubTopics(iSubtopic).topicName);
    classDefinedSubTopics(iSubtopic).subtopicName = strtrim(classDefinedSubTopics(iSubtopic).subtopicName);
    topicIndex = classDefinedTopicsIndex(lower(classDefinedSubTopics(iSubtopic).topicName));
    classDefinedTopics(topicIndex).subtopics(end+1).subtopicName = classDefinedSubTopics(iSubtopic).subtopicName;
    classDefinedTopics(topicIndex).subtopicsIndex(lower(classDefinedSubTopics(iSubtopic).subtopicName)) = length(classDefinedTopics(topicIndex).subtopics);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Step 2) take the method/property/dimension/etc metadata sort those into
% the appropriate topics.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

methodNames = keys(metadataNameMap);
for i=1:length(methodNames)
    metadata = metadataNameMap(methodNames{i});

    % Force everything to have a subtopic
    if ~isfield(metadata,'topic') || isempty(metadata.topic)
        metadata.topic = 'Other';
    end

    if ~isfield(metadata,'subtopic') || isempty(metadata.subtopic)
        metadata.subtopic = 'Other';
    end

    % If this topic isn't created, create it!
    if ~isKey(classDefinedTopicsIndex,lower(metadata.topic))
        classDefinedTopics(end+1).topicName = metadata.topic;
        topicIndex = length(classDefinedTopics);
        classDefinedTopicsIndex(lower(metadata.topic)) = topicIndex;
        classDefinedTopics(topicIndex).subtopics = struct('subtopicName',[],'methodNames',{});
        classDefinedTopics(topicIndex).subtopicsIndex = containers.Map;
    end
    topicIndex = classDefinedTopicsIndex(lower(metadata.topic));

    % If the subtopic isn't created, create it!
    if ~isKey(classDefinedTopics(topicIndex).subtopicsIndex,lower(metadata.subtopic))
        classDefinedTopics(topicIndex).subtopics(end+1).subtopicName = metadata.subtopic;
        classDefinedTopics(topicIndex).subtopicsIndex(lower(metadata.subtopic)) = length(classDefinedTopics(topicIndex).subtopics);
    end
    subtopicIndex = classDefinedTopics(topicIndex).subtopicsIndex(lower(metadata.subtopic));
    
    % And finally, add this method name to the appropriate subtopic index
    classDefinedTopics(topicIndex).subtopics(subtopicIndex).methodNames{end+1} = metadata.name;
end