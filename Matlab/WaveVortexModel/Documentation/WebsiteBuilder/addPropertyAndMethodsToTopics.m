function addPropertyAndMethodsToTopics(rootTopic,metadataNameMap)
% Adds methods and property metadata to the topic list
%
% - Parameter rootTopic: the rootTopic returned by topicsFromClass
% - Parameter metadataNameMap: map with keys of method names and values of normalized metadata structs. 
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
    if isempty(metadata.topic)
        otherTopic.addMethod(metadata);
        continue;
    end

    % if the topic does not exist, create it
    if isempty(rootTopic.subtopicWithName(metadata.topic))
        rootTopic.addSubtopic(Topic(metadata.topic));
    end
    topic = rootTopic.subtopicWithName(metadata.topic);

    % if there is no subtopic, assign to the topic, and then exit
    if isempty(metadata.subtopic)
        topic.addMethod(metadata);
        continue;
    end

    % if the subtopic does not exist, create it
    if isempty(topic.subtopicWithName(metadata.subtopic))
        topic.addSubtopic(Topic(metadata.subtopic));
    end
    subtopic = topic.subtopicWithName(metadata.subtopic);

    % if there is no subsubtopic, assign to the subtopic, and then exit
    if isempty(metadata.subsubtopic)
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