function metadata = ExtractMethodMetadata(mp)
% Extract documentation from method or property (mp) metadata.
%
% Builds a struct that may have the following keys:
% - topic
% - subtopic
% - declaration
% - shortDescription
% - detailedDescription
% - parameters
% - returns
% - className
% - name

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

metadata = ExtractMetadataFromDetailedDescription(mp.DetailedDescription);
metadata.name = mp.Name;
metadata.className = mp.DefiningClass.Name;
metadata.shortDescription = mp.Description;

end