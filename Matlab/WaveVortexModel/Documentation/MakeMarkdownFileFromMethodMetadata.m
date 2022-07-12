function MakeMarkdownFileFromMethodMetadata(path,metadata,pageNumber)
fileID = fopen(path,'w');

fprintf(fileID,'---\nlayout: default\ntitle: %s\nparent: %s\ngrand_parent: Classes\nnav_order: %d\n---\n\n',metadata.name,metadata.className,pageNumber);

fprintf(fileID,'#  %s\n',metadata.name);
fprintf(fileID,'\n%s\n',metadata.shortDescription);
fprintf(fileID,'\n\n---\n\n');
if isfield(metadata,'declaration') && ~isempty(metadata.declaration)
    fprintf(fileID,'## Declaration\n');
    fprintf(fileID,'```matlab\n%s\n```\n',metadata.declaration);
end

if isfield(metadata,'parameters') && ~isempty(metadata.parameters)
    fprintf(fileID,'## Parameters\n');
    for iParameter=1:length(metadata.parameters)
        fprintf(fileID,'+ `%s` %s\n',metadata.parameters(iParameter).name,metadata.parameters(iParameter).description);
    end
    fprintf(fileID,'\n');
end

if isfield(metadata,'detailedDescription') && ~isempty(metadata.detailedDescription)
    fprintf(fileID,'## Discussion\n%s\n',metadata.detailedDescription);
end

fclose(fileID);
end