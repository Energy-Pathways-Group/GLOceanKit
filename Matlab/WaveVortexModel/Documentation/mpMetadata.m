function metadata = mpMetadata(mp)
metadata = [];

% First check if we even want to create documentation for this particular
% property or method.
if isa(mp,'meta.method')
    if strcmp(mp.DefiningClass.Name,'handle') || ~strcmp(mp.Access,'public') || (mp.Hidden == true)
        return;
    end
elseif isa(mp,'meta.property')
    if strcmp(mp.DefiningClass.Name,'handle') || ~strcmp(mp.GetAccess,'public')
        return;
    end
end

% Check out https://regexr.com for testing these regex.
topicExpression = '- topic:(?<topic>[^\r\n]+)(?:$|\n)';
declarationExpression = '- declaration:(?<declaration>[^\r\n]+)(?:$|\n)';
parameterExpression = '- parameter (?<name>[^:]+):(?<description>[^\r\n]+)(?:$|\n)';
leadingWhitespaceExpression = '^[ \t]+';

% Capture the topic annotation, then remove it
detailedDescription = mp.DetailedDescription;
matchStr = regexpi(detailedDescription,topicExpression,'names');
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