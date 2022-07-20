function ClassDocGenerator(className,classDocumentationFolder)
% Generates documentation for a class from the class metadata
%
% Specifically, this script takes the class metadata and creates a markdown
% page for that class using the class's detailed description as the primary
% text, followed by an index, sorted by topic, to all the properties and
% methods of the class.
%
% The metadata for the property and methods is extracted from the property
% and method metadata (from Matlab's metaclass) as well. Each property and
% method gets its own page.

targetFolder = sprintf('%s/%s',classDocumentationFolder,lower(className));
mc = meta.class.fromName(className);

% extract topics and the detailed description (minus those topics)
topicExpression = '- topic:([ \t]*)(?<name>[^\r\n]+)(?:$|\n)';
topics = regexpi(mc.DetailedDescription,topicExpression,'names');
classDetailedDescription = regexprep(mc.DetailedDescription,topicExpression,'','ignorecase');

% Capture metadata from all the public methods and properties
methodAndPropertiesByTopic = ExtractMethodMetadataByTopicFromMetaClass(mc);

% Now create a consolidated list of topics (stored as key names)
mpkeys = methodAndPropertiesByTopic.keys;
commonKeys = intersect({topics(:).name},mpkeys,'stable');
extraKeys = setdiff(mpkeys,{topics(:).name});
allKeys = union(commonKeys,extraKeys,'stable');

% Make a folder for all the class contents
if ~exist(targetFolder,'dir')
    mkdir(targetFolder);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Make the index/table of contents
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileID = fopen(sprintf('%s/index.md',targetFolder),'w');
fprintf(fileID,'---\nlayout: default\ntitle: %s\nparent: Classes\nhas_children: false\nhas_toc: false\n---\n\n',className);
fprintf(fileID,'#  %s\n',className);
fprintf(fileID,'\n%s\n\n',mc.Description);
if ~isempty(mc.DetailedDescription)
    fprintf(fileID,'## Overview\n%s\n',classDetailedDescription);
end
fprintf(fileID,'\n\n## Topics\n');

for iKey=1:length(allKeys)
    topicName = allKeys{iKey};
    mdArray = methodAndPropertiesByTopic(topicName);
    %     fprintf(fileID,'+ [%s](#%s)\n',topicName,replace(lower(topicName),' ','-'));
    fprintf(fileID,'+ %s\n',topicName);
    for i=1:length(mdArray)
        fprintf(fileID,'  + [`%s`](/classes/%s/%s.html) ',mdArray(i).name,lower(className),lower(mdArray(i).name));
        fprintf(fileID,'%s\n',mdArray(i).shortDescription);
    end
end
fprintf(fileID,'\n\n---');

fclose(fileID);

iPageNumber = 0;
for iKey=1:length(allKeys)
    topicName = allKeys{iKey};
    mdArray = methodAndPropertiesByTopic(topicName);
    for i=1:length(mdArray)
        iPageNumber = iPageNumber+1;
        path = sprintf('%s/%s.md',targetFolder,lower(mdArray(i).name));
        MakeMarkdownFileFromMethodMetadata(path,mdArray(i),iPageNumber);
    end
end
end