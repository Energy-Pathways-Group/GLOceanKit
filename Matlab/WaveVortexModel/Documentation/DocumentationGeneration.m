className = 'WaveVortexModel';
mc = meta.class.fromName(className);


methodAndPropertiesByTopic = containers.Map;

% Check out https://regexr.com for testing these regex.
topicExpression = '- topic:(?<topic>[^\r\n]+)(?:$|\n)';
declarationExpression = '- declaration:(?<declaration>[^\r\n]+)(?:$|\n)';
parameterExpression = '- parameter (?<name>[^:]+):(?<description>[^\r\n]+)(?:$|\n)';
leadingWhitespaceExpression = '^[ \t]+';
for i=1:length(mc.MethodList)
    if strcmp(mc.MethodList(i).DefiningClass.Name,'handle') || ~strcmp(mc.MethodList(i).Access,'public') || (mc.MethodList(i).Hidden == true)
        continue;
    end

    mthd = mc.MethodList(i);

    % Capture the topic annotation, then remove it
    detailedDescription = mthd.DetailedDescription;
    matchStr = regexpi(detailedDescription,topicExpression,'names');
    detailedDescription = regexprep(detailedDescription,topicExpression,'','ignorecase');

    if ~isempty(matchStr)
        topicName = regexprep(matchStr.topic,leadingWhitespaceExpression,'');
    else
        topicName = 'Other';
    end

    if isKey(methodAndPropertiesByTopic,topicName)
        methods = methodAndPropertiesByTopic(topicName);
    else
        methods = [];
    end

    methods(end+1).name = mthd.Name;

    % Capture all parameters, then remove the annotations
    matchStr = regexpi(detailedDescription,parameterExpression,'names');
    detailedDescription = regexprep(detailedDescription,parameterExpression,'','ignorecase');
    if ~isempty(matchStr)
        methods(end).parameters = matchStr;
    end

    % Capture any declarations made, then remove the annotation
    matchStr = regexpi(detailedDescription,declarationExpression,'names');
    detailedDescription = regexprep(detailedDescription,declarationExpression,'','ignorecase');
    if ~isempty(matchStr)
        methods(end).declaration = matchStr.declaration;
    end

    methods(end).shortDescription = mthd.Description;
    methods(end).detailedDescription = regexprep(detailedDescription,leadingWhitespaceExpression,'');

    methodAndPropertiesByTopic(topicName) = methods;

%     fprintf('%d: %s\n',i,mc.MethodList(i).Name)
end

for i=1:length(mc.PropertyList)
    if strcmp(mc.PropertyList(i).DefiningClass.Name,'handle') || ~strcmp(mc.PropertyList(i).GetAccess,'public')
        continue;
    end

    mthd = mc.PropertyList(i);

    % Capture the topic annotation, then remove it
    detailedDescription = mthd.DetailedDescription;
    matchStr = regexpi(detailedDescription,topicExpression,'names');
    detailedDescription = regexprep(detailedDescription,topicExpression,'','ignorecase');

    if ~isempty(matchStr)
        topicName = regexprep(matchStr.topic,leadingWhitespaceExpression,'');
    else
        topicName = 'Other';
    end

    if isKey(methodAndPropertiesByTopic,topicName)
        methods = methodAndPropertiesByTopic(topicName);
    else
        methods = [];
    end

    methods(end+1).name = mthd.Name;

    % Capture all parameters, then remove the annotations
    matchStr = regexpi(detailedDescription,parameterExpression,'names');
    detailedDescription = regexprep(detailedDescription,parameterExpression,'','ignorecase');
    if ~isempty(matchStr)
        methods(end).parameters = matchStr;
    end

    % Capture any declarations made, then remove the annotation
    methods(end).declaration = regexpi(detailedDescription,declarationExpression,'match');
    detailedDescription = regexprep(detailedDescription,declarationExpression,'','ignorecase');

    methods(end).shortDescription = mthd.Description;
    methods(end).detailedDescription = regexprep(detailedDescription,leadingWhitespaceExpression,'');

    methodAndPropertiesByTopic(topicName) = methods;
end

% Make a folder for all the class contents
if ~exist(sprintf('%s',lower(className)),'dir')
    mkdir(sprintf('%s',lower(className)))
end

%%%%%%%%%%%%%%%%%%%%%
%
% Make the index/table of contents
%
fileID = fopen(sprintf('%s/index.md',lower(className)),'w');
fprintf(fileID,'---\nlayout: default\ntitle: %s\nparent: Classes\nhas_children: true\n---\n\n',className);
fprintf(fileID,'#  %s\n',className);
fprintf(fileID,'\n%s\n\n',mc.Description);
if ~isempty(mc.DetailedDescription)
    fprintf(fileID,'## Discussion\n%s\n',mc.DetailedDescription);
end
fprintf(fileID,'\n\n## Topics\n');
mpkeys = methodAndPropertiesByTopic.keys;
for iKey=1:length(mpkeys)
    topicName = mpkeys{iKey};
    methods = methodAndPropertiesByTopic(topicName);
    fprintf(fileID,'+ [%s](#%s)\n',topicName,replace(lower(topicName),' ','-'));
    for i=1:length(methods)
        fprintf(fileID,'  + [`%s`](/classes/%s/%s.html) ',methods(i).name,lower(className),lower(methods(i).name));
        fprintf(fileID,'%s\n',methods(i).shortDescription);
    end
end
fprintf(fileID,'\n\n---');

fclose(fileID);


for iKey=1:length(mpkeys)
    topicName = mpkeys{iKey};
    methods = methodAndPropertiesByTopic(topicName);
    for i=1:length(methods)
        fileID = fopen(sprintf('%s/%s.md',lower(className),lower(methods(i).name)),'w');
        
        fprintf(fileID,'---\nlayout: default\ntitle: %s\nparent: %s\ngrand_parent: Classes\n---\n\n',methods(i).name,className);

        fprintf(fileID,'#  %s\n',methods(i).name);
        fprintf(fileID,'\n%s\n',methods(i).shortDescription);
        fprintf(fileID,'\n\n---\n\n');
        if isfield(methods(i),'declaration') && ~isempty(methods(i).declaration)
            fprintf(fileID,'## Declaration\n');
            fprintf(fileID,'```matlab\n%s\n```\n',methods(i).declaration);
        end

        if isfield(methods(i),'parameters') && ~isempty(methods(i).parameters)
            fprintf(fileID,'## Parameters\n');
            for iParameter=1:length(methods(i).parameters)
                fprintf(fileID,'+ `%s` %s\n',methods(i).parameters(iParameter).name,methods(i).parameters(iParameter).description);
            end
            fprintf(fileID,'\n');
        end

        if isfield(methods(i),'detailedDescription') && ~isempty(methods(i).detailedDescription)
            fprintf(fileID,'## Discussion\n%s\n',methods(i).detailedDescription);
        end

        fclose(fileID);
    end
end


%%%%%%%%%%%%%%%%%%%%%
%
% Now write out the help for each method
%
% for iKey=1:length(mpkeys)
%     topicName = mpkeys{iKey};
%     fprintf(fileID,'\n\n## %s\n\n---',topicName);
%     methods = methodAndPropertiesByTopic(topicName);
%     for i=1:length(methods)
%         fprintf(fileID,'\n\n### `%s`\n',methods(i).name);
%         fprintf(fileID,'\n%s\n',methods.shortDescription);
%         if isfield(methods,'declaration') && ~isempty(methods.declaration)
%             
%         end
%         fprintf(fileID,'  + [`%s`](#%s)\n',methods(i).name,lower(methods(i).name));
%     end
% end



%     str = regexprep(str,expression,'','ignorecase');