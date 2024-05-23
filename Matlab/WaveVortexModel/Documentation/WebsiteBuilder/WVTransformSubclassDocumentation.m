classdef WVTransformSubclassDocumentation < ClassDocumentation
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    methods

        function self = WVTransformSubclassDocumentation(name,options)
            arguments
                name
                options.buildFolder % the folder where we are dumping everything on the local hard drive. This will become the *root* website folder
                options.websiteFolder % the folder relative to the root website folder
                options.parent = []
                options.grandparent = []
                options.nav_order = []
                options.excludedSuperclasses = {'handle'};
            end
            superClassOptions = namedargs2cell(options);
            self@ClassDocumentation(name,superClassOptions{:});

            classMetadata = meta.class.fromName(self.name);
            for iClass=1:length(classMetadata.SuperclassList)
                switch classMetadata.SuperclassList(iClass).Name
                    case 'WVStratifiedFlow'
                        props = WVStratifiedFlow.propertyAnnotationsForStratifiedFlow();
                        for iDim=1:length(props)
                            metadata = MethodDocumentation(props(iDim).name);
                            metadata.addMetadataFromDetailedDescription(props(iDim).detailedDescription);
                            metadata.definingClassName = self.name;
                            metadata.addDeclaringClass(self.name);
                            metadata.shortDescription = props(iDim).description;
                            metadata.units = props(iDim).units;
                            metadata.dimensions = props(iDim).dimensions;
                            metadata.isComplex = props(iDim).isComplex;
                            metadata.functionType = FunctionType.transformProperty;
                            self.addMethodDocumentation(metadata);
                        end

                        methods = WVStratifiedFlow.methodAnnotationsForStratifiedFlow;
                        for iMethod=1:length(methods)
                            metadata = MethodDocumentation(methods(iMethod).name);
                            metadata.addMetadataFromDetailedDescription(methods(iMethod).detailedDescription);
                            metadata.shortDescription = methods(iMethod).description;
                            metadata.definingClassName = self.name;
                            metadata.addDeclaringClass(self.name);
                            metadata.functionType = FunctionType.instanceMethod;
                            self.addMethodDocumentation(metadata);
                        end
                end
            end

        end
    end

    methods (Static)
        function classDocumentation = classDocumentationFromClassNames(nameArray,options)
            % initialize several ClassDocumentations at the same time. This
            % is useful when the documentation is at the same level
            % (sharing the same parent and grandparent). The documentation
            % will be assigned a nav_order.
            arguments
                nameArray
                options.buildFolder % the folder where we are dumping everything on the local hard drive. This will become the *root* website folder
                options.websiteFolder % the folder relative to the root website folder
                options.parent = []
                options.grandparent = []
                options.excludedSuperclasses = {'handle'};
            end

            classDocumentation = WVTransformSubclassDocumentation.empty(length(nameArray),0);
            for iName=1:length(nameArray)
                classDocumentation(iName) = WVTransformSubclassDocumentation(nameArray{iName},buildFolder=options.buildFolder,websiteFolder=options.websiteFolder,parent=options.parent,grandparent=options.grandparent,nav_order=iName,excludedSuperclasses = options.excludedSuperclasses);
            end
        end
    end
end