function metadataNameMap = metadataFromClassPropertiesAndMethods(mc)
% Capture metadata from all the public methods and properties of a class
%
% This function ultimately calls -ExtractMetadataFromDetailedDescription,
% but first sorts through the available methods and properties to find the
% right ones.
%
% Builds a struct that may have the following keys:
% - topic
% - subtopic
% - subsubtopic
% - declaration
% - shortDescription
% - detailedDescription
% - parameters
% - returns
% - className
% - name
%
% - Parameter mc: the detailed description 
% - Returns metadataNameMap: containers.Map with method names as keys and metadata structures as values.
arguments
    mc meta.class
end

metadataNameMap = containers.Map;
for i=1:length(mc.MethodList)
    metadata = ExtractMethodMetadata(mc.MethodList(i),mc.Name);
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
    metadata = ExtractMethodMetadata(mc.PropertyList(i),mc.Name);
    if ~isempty(metadata)
        metadata.functionType = FunctionType.instanceProperty;
        metadataNameMap(metadata.name) = metadata;
    end
end
end

function metadata = ExtractMethodMetadata(mp,className)
% Extract documentation from method or property (mp) metadata.
metadata = [];

% Don't create documentation if this is a method defined in the superclass
% This initial check does not work, because if the subclass re-defines a
% method, then it counts the subclass as the "DefiningClass". But, for
% documentation purposes, we really don't want that.
if ~strcmp(mp.DefiningClass.Name,className)
    return;
end
if isMethodDefinedInSuperclass(mp)
    return;
end

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

metadata = MethodAnnotation(mp.Name);
metadata.addMetadataFromDetailedDescription(mp.DetailedDescription);
metadata.className = mp.DefiningClass.Name;
metadata.shortDescription = mp.Description;

end

function bool = isMethodDefinedInSuperclass(mp)
bool = 0;
if isempty(mp.DefiningClass.SuperclassList)
    return;
end
for i=1:length(mp.DefiningClass.SuperclassList(1).MethodList)
    if strcmp(mp.DefiningClass.SuperclassList(1).MethodList(i).Name,mp.Name)
        bool = 1;
        return;
    end
end
end