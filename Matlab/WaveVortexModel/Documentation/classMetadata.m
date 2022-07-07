function topics = classMetadata(mc)

% Check out https://regexr.com for testing these regex.
topicExpression = '- topic:(?<topic>[^\r\n]+)(?:$|\n)';
topics = regexpi(mc.DetailedDescription,topicExpression,'names');

methodAndPropertiesByTopic = containers.Map;



% Capture the topic annotation, then remove it
detailedDescription = mp.DetailedDescription;

detailedDescription = regexprep(detailedDescription,topicExpression,'','ignorecase');

if ~isempty(matchStr)
    topicName = regexprep(matchStr.topic,leadingWhitespaceExpression,'');
else
    topicName = 'Other';
end

metadata.topic = topicName;
metadata.name = mp.Name;

% Capture all parameters, then remove the annotations
metadata.parameters = regexpi(detailedDescription,parameterExpression,'names');
detailedDescription = regexprep(detailedDescription,parameterExpression,'','ignorecase');

% Capture any declarations made, then remove the annotation
matchStr = regexpi(detailedDescription,declarationExpression,'names');
detailedDescription = regexprep(detailedDescription,declarationExpression,'','ignorecase');
if ~isempty(matchStr)
    metadata.declaration = matchStr.declaration;
else
    metadata.declaration = [];
end

metadata.shortDescription = mp.Description;
metadata.detailedDescription = regexprep(detailedDescription,leadingWhitespaceExpression,'');

end