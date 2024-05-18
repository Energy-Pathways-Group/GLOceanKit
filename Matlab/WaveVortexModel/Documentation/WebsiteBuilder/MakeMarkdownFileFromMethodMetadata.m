function MakeMarkdownFileFromMethodMetadata(path,metadata,pageNumber)
fileID = fopen(path,'w');

fprintf(fileID,'---\nlayout: default\ntitle: %s\nparent: %s\ngrand_parent: Classes\nnav_order: %d\nmathjax: true\n---\n\n',metadata.name,metadata.className,pageNumber);

fprintf(fileID,'#  %s\n',metadata.name);
fprintf(fileID,'\n%s\n',metadata.shortDescription);

fprintf(fileID,'\n\n---\n\n');

if (metadata.functionType == FunctionType.transformProperty || metadata.functionType == FunctionType.stateVariable)
    fprintf(fileID,'## Description\n');
    if metadata.isComplex == 1
        str = 'Complex valued ';
    else
        str = 'Real valued ';
    end

    if metadata.functionType == FunctionType.transformProperty
        str = strcat(str,' transform property ');
    else
        str = strcat(str,' state variable ');
    end

    if isempty(metadata.dimensions)
        str = strcat(str,' with no dimensions and ');
    elseif length(metadata.dimensions) == 1 
        str = strcat(str,' with dimension %s and ',metadata.dimensions{1});
    else
        str = strcat(str,' with dimensions $$(');
        for iDim=1:(length(metadata.dimensions)-1)
            str = strcat(str,metadata.dimensions{iDim},',');
        end
        str = strcat(str,metadata.dimensions{end},')$$ and');
    end

    if isempty(metadata.units)
        str = strcat(str,' no units.\n\n');
    else
        str = strcat(str,' units of $$',metadata.units,'$$.\n\n');
    end
    fprintf(fileID,str);
end

% ## Description
% (Real/Complex) valued (transform property/state variable) with (no
% dimensions/dimension {x,y,z}) and (no units/units of xxx).



if ~isempty(metadata.declaration)
    fprintf(fileID,'## Declaration\n');
    fprintf(fileID,'```matlab\n%s\n```\n',metadata.declaration);
end

if ~isempty(metadata.parameters)
    fprintf(fileID,'## Parameters\n');
    for iParameter=1:length(metadata.parameters)
        fprintf(fileID,'+ `%s` %s\n',metadata.parameters(iParameter).name,metadata.parameters(iParameter).description);
    end
    fprintf(fileID,'\n');
end

if ~isempty(metadata.returns)
    fprintf(fileID,'## Returns\n');
    for iReturn=1:length(metadata.returns)
        fprintf(fileID,'+ `%s` %s\n',metadata.returns(iReturn).name,metadata.returns(iReturn).description);
    end
    fprintf(fileID,'\n');
end

if ~isempty(metadata.detailedDescription)
    fprintf(fileID,'## Discussion\n%s\n',metadata.detailedDescription);
end

fclose(fileID);
end