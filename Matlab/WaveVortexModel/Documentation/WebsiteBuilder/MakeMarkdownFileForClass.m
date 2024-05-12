function MakeMarkdownFileForClass(options)
arguments
    options.path
    options.websiteFolder
    options.className
    options.classDetailedDescription
    options.classDefinedTopics
    options.metadataNameMap
    options.parent = []
    options.grandparent = []
    options.nav_order = []
end
className = options.className;
classDefinedTopics = options.classDefinedTopics;
metadataNameMap = options.metadataNameMap;

mc = meta.class.fromName(className);

declarationExpression = '- declaration:(?<declaration>[^\r\n]+)(?:$|\n)';
matchStr = regexpi(options.classDetailedDescription,declarationExpression,'names');
classDetailedDescription = regexprep(options.classDetailedDescription,declarationExpression,'','ignorecase');
if ~isempty(matchStr)
    declaration = matchStr.declaration;
else
    declaration = [];
end

fileID = fopen(options.path,'w');
fprintf(fileID,'---\nlayout: default\ntitle: %s\nhas_children: false\nhas_toc: false\nmathjax: true\n',className);
if ~isempty(options.parent)
    fprintf(fileID,'parent: %s\n',options.parent);
end
if ~isempty(options.grandparent)
    fprintf(fileID,'grand_parent: %s\n',options.grandparent);
end
if ~isempty(options.nav_order)
    fprintf(fileID,'nav_order: %d\n',options.nav_order);
end
fprintf(fileID,'---\n\n');

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

% do not call the method directly, because we do not want to print the Root
% category name
for iSubtopic = 1:length(classDefinedTopics.subtopics)
    writeMarkdownForTopic(classDefinedTopics.subtopics(iSubtopic),fileID,metadataNameMap,'',options.websiteFolder);
end

fprintf(fileID,'\n\n---');

fclose(fileID);
end

function writeMarkdownForTopic(topic,fileID,metadataNameMap,indentLevel,websiteFolder)

if isempty(topic.methodNames) && isempty(topic.subtopics)
    return
end
fprintf(fileID,'%s+ %s\n',indentLevel,topic.name);
for methodIndex = 1:length(topic.methodNames)
    fprintf(fileID,'%s  + [`%s`](/%s/%s.html) ',indentLevel,topic.methodNames{methodIndex},websiteFolder,lower(topic.methodNames{methodIndex}));
    fprintf(fileID,'%s\n',metadataNameMap(topic.methodNames{methodIndex}).shortDescription);
end
for iSubtopic = 1:length(topic.subtopics)
    writeMarkdownForTopic(topic.subtopics(iSubtopic),fileID,metadataNameMap,[indentLevel,'  '],websiteFolder);
end

end