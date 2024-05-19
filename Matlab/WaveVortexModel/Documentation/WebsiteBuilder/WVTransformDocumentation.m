classdef WVTransformDocumentation < ClassDocumentation
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    methods

        function initializeMethodDocumentation(self)
            initializeMethodDocumentation@ClassDocumentation(self);

            % Overwrite the automatically extracted metadata, with our version
            dims = WVTransform.defaultDimensionAnnotations();
            for iDim=1:length(dims)
                metadata = MethodDocumentation(dims(iDim).name);
                metadata.addMetadataFromDetailedDescription(dims(iDim).detailedDescription);
                metadata.className = self.name;
                metadata.shortDescription = dims(iDim).description;
                metadata.units = dims(iDim).units;
                metadata.isComplex = 0;
                metadata.functionType = FunctionType.transformDimension;
                self.addMethodDocumentation(metadata);
            end

            % Overwrite the automatically extracted metadata, with our version
            props = WVTransform.defaultPropertyAnnotations();
            for iDim=1:length(props)
                metadata = MethodDocumentation(props(iDim).name);
                metadata.addMetadataFromDetailedDescription(props(iDim).detailedDescription);
                metadata.className = self.name;
                metadata.shortDescription = props(iDim).description;
                metadata.units = props(iDim).units;
                metadata.dimensions = props(iDim).dimensions;
                metadata.isComplex = props(iDim).isComplex;
                metadata.functionType = FunctionType.transformProperty;
                self.addMethodDocumentation(metadata);
            end

            % Overwrite the automatically extracted metadata, with our version
            ops = WVTransform.defaultOperations();
            for iOp=1:length(ops)
                for iVar=1:length(ops(iOp).outputVariables)
                    stateVar = ops(iOp).outputVariables(iVar);

                    metadata = MethodDocumentation(stateVar.name);
                    metadata.addMetadataFromDetailedDescription(stateVar.detailedDescription);
                    metadata.className = self.name;
                    metadata.shortDescription = stateVar.description;
                    metadata.units = stateVar.units;
                    metadata.dimensions = stateVar.dimensions;
                    metadata.isComplex = stateVar.isComplex;
                    metadata.functionType = FunctionType.stateVariable;
                    self.addMethodDocumentation(metadata);
                end
            end

            methods = WVTransform.defaultMethodAnnotations;
            for iMethod=1:length(methods)
                metadata = MethodDocumentation(methods(iMethod).name);
                metadata.addMetadataFromDetailedDescription(methods(iMethod).detailedDescription);
                metadata.shortDescription = methods(iMethod).description;
                metadata.className = self.name;
                metadata.functionType = FunctionType.instanceMethod;
                self.addMethodDocumentation(metadata);
            end
        end
    end
end