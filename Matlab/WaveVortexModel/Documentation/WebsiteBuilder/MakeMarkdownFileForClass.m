function MakeMarkdownFileForClass(path,className,classDetailedDescription,classDefinedTopics, metadataNameMap,parentName,parentFolder)

mc = meta.class.fromName(className);

declarationExpression = '- declaration:(?<declaration>[^\r\n]+)(?:$|\n)';
matchStr = regexpi(classDetailedDescription,declarationExpression,'names');
classDetailedDescription = regexprep(classDetailedDescription,declarationExpression,'','ignorecase');
if ~isempty(matchStr)
    declaration = matchStr.declaration;
else
    declaration = [];
end

fileID = fopen(path,'w');
fprintf(fileID,'---\nlayout: default\ntitle: %s\nparent: %s\nhas_children: false\nhas_toc: false\nmathjax: true\n---\n\n',className,parentName);
fprintf(fileID,'#  %s\n',className);
fprintf(fileID,'\n%s\n',mc.Description);

fprintf(fileID,'\n\n---\n\n');

if ~isempty(declaration)
    fprintf(fileID,'## Declaration\n');
    %fprintf(fileID,'```matlab\n%s\n```\n\n',declaration);
    % unfortunately we cannot directly put links into a markdown code block, so
    % we have to manually extract the url ourselves, and then write our own
    % html. so annoying!!
    mdurlExpression = '\[(?<name>[^\[]+)\]\((?<url>.*)\)';
    matchStr = regexpi(declaration,mdurlExpression,'names');
    if ~isempty(matchStr)
        linkString = strcat('<a href="',matchStr.url,'" title="',matchStr.name,'">',matchStr.name,'</a>');
        declaration = regexprep(declaration,mdurlExpression,linkString,'ignorecase');
    end
    prestring = '<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>';
    poststring = '</code></pre></div></div>';
    fprintf(fileID,'\n%s\n\n',strcat(prestring,strip(declaration),poststring));
end

if ~isempty(mc.DetailedDescription)
    fprintf(fileID,'## Overview\n%s\n',classDetailedDescription);
end
fprintf(fileID,'\n\n## Topics\n');

for topicIndex = 1:length(classDefinedTopics)
    % If the 'Other' subtopic isn't there, it means no methods actually
    % used this topic.
%     if ~isKey(classDefinedTopics(topicIndex).subtopicsIndex,lower('Other'))
%         continue;
%     end

    fprintf(fileID,'+ %s\n',classDefinedTopics(topicIndex).topicName);

    % First do the 'other' subtopics immediately below
    if isKey(classDefinedTopics(topicIndex).subtopicsIndex,lower('Other'))
        otherSubtopicIndex = classDefinedTopics(topicIndex).subtopicsIndex(lower('Other'));

        subtopic = classDefinedTopics(topicIndex).subtopics(otherSubtopicIndex);
        for methodIndex = 1:length(subtopic.methodNames)
            fprintf(fileID,'  + [`%s`](/%s/%s/%s.html) ',subtopic.methodNames{methodIndex},parentFolder,lower(className),lower(subtopic.methodNames{methodIndex}));
            fprintf(fileID,'%s\n',metadataNameMap(subtopic.methodNames{methodIndex}).shortDescription);
        end
    else
        otherSubtopicIndex = 0;
    end
    for subtopicIndex = 1:length(classDefinedTopics(topicIndex).subtopics)
        if subtopicIndex == otherSubtopicIndex
            continue;
        end
        subtopic = classDefinedTopics(topicIndex).subtopics(subtopicIndex);
        fprintf(fileID,'  + %s\n',subtopic.subtopicName);
        for methodIndex = 1:length(subtopic.methodNames)
            fprintf(fileID,'    + [`%s`](/%s/%s/%s.html) ',subtopic.methodNames{methodIndex},parentFolder,lower(className),lower(subtopic.methodNames{methodIndex}));
            fprintf(fileID,'%s\n',metadataNameMap(subtopic.methodNames{methodIndex}).shortDescription);
        end
    end
end

fprintf(fileID,'\n\n---');

fclose(fileID);
end