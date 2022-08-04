function metadataNameMap = ExtractMetadataFromClassPropertiesAndMethods(mc)
% Capture metadata from all the public methods and properties of a class
%
% - Parameter mc: the detailed description 
% - Returns metadataNameMap: containers.Map with method names as keys and metadata structures as values.
arguments
    mc meta.class
end

metadataNameMap = containers.Map;
for i=1:length(mc.MethodList)
    metadata = ExtractMethodMetadata(mc.MethodList(i));
    if ~isempty(metadata)
        if mc.MethodList(i).Static == 1
            metadata.functionType = FunctionType.staticMethod;
        elseif mc.MethodList(i).Abstract == 1
            metadata.functionType = FunctionType.abstractMethod;
        else
            metadata.functionType = FunctionType.instanceMethod;
        end
        metadataNameMap(metadata.name) = metadata;
    end
end
for i=1:length(mc.PropertyList)
    metadata = ExtractMethodMetadata(mc.PropertyList(i));
    if ~isempty(metadata)
        metadata.functionType = FunctionType.instanceProperty;
        metadataNameMap(metadata.name) = metadata;
    end
end
end

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