function ClassDocGenerator(className,classDocumentationFolder)
targetFolder = sprintf('%s/%s',classDocumentationFolder,lower(className));
mc = meta.class.fromName(className);

topicExpression = '- topic:([ \t]*)(?<name>[^\r\n]+)(?:$|\n)';
topics = regexpi(mc.DetailedDescription,topicExpression,'names');
classDetailedDescription = regexprep(mc.DetailedDescription,topicExpression,'','ignorecase');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Capture metadata from all the public methods and properties
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

methodAndPropertiesByTopic = containers.Map;
for i=1:length(mc.MethodList)
    metadata = mpMetadata(mc.MethodList(i));
    if ~isempty(metadata)
        if isKey(methodAndPropertiesByTopic,metadata.topic)
            mdArray = methodAndPropertiesByTopic(metadata.topic);
            mdArray(end+1) = metadata;
        else
            mdArray = metadata;
        end
        methodAndPropertiesByTopic(metadata.topic) = mdArray;
    end
end

for i=1:length(mc.PropertyList)
    metadata = mpMetadata(mc.PropertyList(i));
    if ~isempty(metadata)
        if isKey(methodAndPropertiesByTopic,metadata.topic)
            mdArray = methodAndPropertiesByTopic(metadata.topic);
            mdArray(end+1) = metadata;
        else
            mdArray = metadata;
        end
        methodAndPropertiesByTopic(metadata.topic) = mdArray;
    end
end

% Make a folder for all the class contents
if ~exist(targetFolder,'dir')
    mkdir(targetFolder);
end

%%%%%%%%%%%%%%%%%%%%%
%
% Make the index/table of contents
%
fileID = fopen(sprintf('%s/index.md',targetFolder),'w');
fprintf(fileID,'---\nlayout: default\ntitle: %s\nparent: Classes\nhas_children: true\nhas_toc: false\n---\n\n',className);
fprintf(fileID,'#  %s\n',className);
fprintf(fileID,'\n%s\n\n',mc.Description);
if ~isempty(mc.DetailedDescription)
    fprintf(fileID,'## Discussion\n%s\n',classDetailedDescription);
end
fprintf(fileID,'\n\n## Topics\n');
mpkeys = methodAndPropertiesByTopic.keys;

commonKeys = intersect({topics(:).name},mpkeys,'stable');
extraKeys = setdiff(mpkeys,{topics(:).name});
allKeys = union(commonKeys,extraKeys,'stable');
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
        fileID = fopen(sprintf('%s/%s.md',targetFolder,lower(mdArray(i).name)),'w');

        fprintf(fileID,'---\nlayout: default\ntitle: %s\nparent: %s\ngrand_parent: Classes\nnav_order: %d\n---\n\n',mdArray(i).name,className,iPageNumber);

        fprintf(fileID,'#  %s\n',mdArray(i).name);
        fprintf(fileID,'\n%s\n',mdArray(i).shortDescription);
        fprintf(fileID,'\n\n---\n\n');
        if isfield(mdArray(i),'declaration') && ~isempty(mdArray(i).declaration)
            fprintf(fileID,'## Declaration\n');
            fprintf(fileID,'```matlab\n%s\n```\n',mdArray(i).declaration);
        end

        if isfield(mdArray(i),'parameters') && ~isempty(mdArray(i).parameters)
            fprintf(fileID,'## Parameters\n');
            for iParameter=1:length(mdArray(i).parameters)
                fprintf(fileID,'+ `%s` %s\n',mdArray(i).parameters(iParameter).name,mdArray(i).parameters(iParameter).description);
            end
            fprintf(fileID,'\n');
        end

        if isfield(mdArray(i),'detailedDescription') && ~isempty(mdArray(i).detailedDescription)
            fprintf(fileID,'## Discussion\n%s\n',mdArray(i).detailedDescription);
        end

        fclose(fileID);
    end
end
end