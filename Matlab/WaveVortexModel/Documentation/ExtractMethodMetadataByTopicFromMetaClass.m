function mpMap = ExtractMethodMetadataByTopicFromMetaClass(mc)
% Capture metadata from all the public methods and properties
arguments
    mc meta.class
end

mpMap = containers.Map;
for i=1:length(mc.MethodList)
    metadata = ExtractMethodMetadata(mc.MethodList(i));
    if ~isempty(metadata)
        if isKey(mpMap,metadata.topic)
            mdArray = mpMap(metadata.topic);
            mdArray(end+1) = metadata;
        else
            mdArray = metadata;
        end
        mpMap(metadata.topic) = mdArray;
    end
end

for i=1:length(mc.PropertyList)
    metadata = ExtractMethodMetadata(mc.PropertyList(i));
    if ~isempty(metadata)
        if isKey(mpMap,metadata.topic)
            mdArray = mpMap(metadata.topic);
            mdArray(end+1) = metadata;
        else
            mdArray = metadata;
        end
        mpMap(metadata.topic) = mdArray;
    end
end
end