function metadata = ExtractMetadataFromDetailedDescription(detailedDescription)
% Builds a struct that may have the following keys:
% - topic
% - subtopic
% - declaration
% - shortDescription
% - detailedDescription
% - parameters
% - returns

metadata = [];
if isempty(detailedDescription)
    return;
end

% Check out https://regexr.com for testing these regex.
topicExpression = '- topic:([ \t]*)(?<topic>[^\r\n]+)(?:$|\n)';
subtopicExpression = '- topic:([ \t]*)(?<topic>[^—\r\n]+)—([ \t]*)(?<subtopic>[^\r\n]+)(?:$|\n)';
declarationExpression = '- declaration:(?<declaration>[^\r\n]+)(?:$|\n)';
parameterExpression = '- parameter (?<name>[^:]+):(?<description>[^\r\n]+)(?:$|\n)';
returnsExpression = '- returns (?<name>[^:]+):(?<description>[^\r\n]+)(?:$|\n)';
leadingWhitespaceExpression = '^[ \t]+';

% Capture the subtopic annotation, then remove it
matchStr = regexpi(detailedDescription,subtopicExpression,'names');
detailedDescription = regexprep(detailedDescription,subtopicExpression,'','ignorecase');
if ~isempty(matchStr)
    metadata.subtopic = strtrim(matchStr.subtopic);
    metadata.topic = strtrim(matchStr.topic);
end

% Capture the topic annotation, then remove it
matchStr = regexpi(detailedDescription,topicExpression,'names');
detailedDescription = regexprep(detailedDescription,topicExpression,'','ignorecase');
if ~isempty(matchStr)
    topicName = strtrim(matchStr.topic);
else
    topicName = 'Other';
end
if ~isfield(metadata,'topic') || isempty(metadata.topic)
    metadata.topic = topicName;
end

% Capture all parameters, then remove the annotations
metadata.parameters = regexpi(detailedDescription,parameterExpression,'names');
detailedDescription = regexprep(detailedDescription,parameterExpression,'','ignorecase');

% Capture all returns, then remove the annotations
metadata.returns = regexpi(detailedDescription,returnsExpression,'names');
detailedDescription = regexprep(detailedDescription,returnsExpression,'','ignorecase');

% Capture any declarations made, then remove the annotation
matchStr = regexpi(detailedDescription,declarationExpression,'names');
detailedDescription = regexprep(detailedDescription,declarationExpression,'','ignorecase');
if ~isempty(matchStr)
    metadata.declaration = matchStr.declaration;
else
    metadata.declaration = [];
end


metadata.detailedDescription = regexprep(detailedDescription,leadingWhitespaceExpression,'');
end