function [detailedDescription,rootTopic] = topicsFromClass(mc)
% Extracts topics and subtopic from a detailedDescription and creates a structure useful creating an ordered topic index
%
% - Parameter mc: the detailed description 
% - Parameter metadataNameMap: map with keys of method names and values of normalized metadata structs. 
% - Returns detailedDescription: the detailed description without the topic/subtopic metadata
% - Returns classDefinedTopics: a structure array with field "topicName" (char), "subtopics" (struct('subtopicName',[],'methodNames',{})), and "subtopicsIndex" (containers.Map with subtopic name mapping to an order index).
arguments
    mc meta.class
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
