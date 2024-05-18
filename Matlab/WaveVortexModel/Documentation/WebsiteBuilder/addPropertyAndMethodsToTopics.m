function addPropertyAndMethodsToTopics(rootTopic,metadataNameMap)
% Extracts topics and subtopic from a detailedDescription and creates a structure useful creating an ordered topic index
%
% - Parameter mc: the detailed description 
% - Parameter metadataNameMap: map with keys of method names and values of normalized metadata structs. 
% - Returns detailedDescription: the detailed description without the topic/subtopic metadata
% - Returns classDefinedTopics: a structure array with field "topicName" (char), "subtopics" (struct('subtopicName',[],'methodNames',{})), and "subtopicsIndex" (containers.Map with subtopic name mapping to an order index).
arguments
    rootTopic Topic
    metadataNameMap containers.Map
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