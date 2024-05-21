classdef ClassDocumentation < handle
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

    properties
        name string
        pathOfClassFolderRelative
        pathOfClassFolderAbsolute
        parent
        grandparent
        nav_order
        rootTopic

        allMethodDocumentation
        methodDocumentationNameMap
        shouldExcludeAllSuperclasses = 1;

        shortDescription
        detailedDescription
        declaration = []
    end

    properties (SetObservable, AbortSet, Access = public)
        excludedSuperclasses = {'handle'};
    end

    methods
        function self = ClassDocumentation(name,options)
            arguments
                name
                options.buildFolder % the folder where we are dumping everything on the local hard drive. This will become the *root* website folder
                options.websiteFolder % the folder relative to the root website folder
                options.parent = []
                options.grandparent = []
                options.nav_order = []
            end
            self.name = name;
            self.parent = options.parent;
            self.grandparent = options.grandparent;
            self.nav_order = options.nav_order;

            % relative path to a folder for all the class contents
            self.pathOfClassFolderRelative = sprintf('%s/%s',options.websiteFolder,lower(self.name));

            % Make a folder for all the class contents
            self.pathOfClassFolderAbsolute = sprintf('%s/',options.buildFolder,self.pathOfClassFolderRelative);


            mc = meta.class.fromName(self.name);
            
            self.shortDescription = mc.Description;
            self.declaration = ClassDocumentation.declarationFromString(mc.DetailedDescription);

            self.detailedDescription = mc.DetailedDescription;
            self.detailedDescription = ClassDocumentation.trimDeclarationFromString(self.detailedDescription);
            self.detailedDescription = Topic.trimTopicsFromString(self.detailedDescription);

            self.allMethodDocumentation = ClassDocumentation.methodDocumentationFromClass(self.name);
            addlistener(self,'excludedSuperclasses','PostSet',@self.excludedSuperclassesDidChange); 
        end

        function excludedSuperclassesDidChange(self,~,~)
            self.shouldExcludeAllSuperclasses = 0;
        end

        function addMethodDocumentation(self,methodDocumentation)
            if isKey(self.methodDocumentationNameMap,methodDocumentation.name)
                methodDocumentation.pathOfOutputFile = self.methodDocumentationNameMap(methodDocumentation.name).pathOfOutputFile;
                methodDocumentation.pathOfFileOnWebsite = self.methodDocumentationNameMap(methodDocumentation.name).pathOfFileOnWebsite;
            elseif ismember(lower(methodDocumentation.name),lower(keys(self.methodDocumentationNameMap))) && ~isKey(self.methodDocumentationNameMap,methodDocumentation.name)
                methodDocumentation.pathOfOutputFile = sprintf('%s/%s_.md',self.pathOfClassFolderAbsolute,lower(methodDocumentation.name));
                methodDocumentation.pathOfFileOnWebsite = sprintf('/%s/%s_.html',self.pathOfClassFolderRelative,lower(methodDocumentation.name));
            else
                methodDocumentation.pathOfOutputFile = sprintf('%s/%s.md',self.pathOfClassFolderAbsolute,lower(methodDocumentation.name));
                methodDocumentation.pathOfFileOnWebsite = sprintf('/%s/%s.html',self.pathOfClassFolderRelative,lower(methodDocumentation.name));
            end
            self.methodDocumentationNameMap(methodDocumentation.name) = methodDocumentation;
        end

        function initializeMethodDocumentation(self)
            % Capture metadata from all the public methods and properties of a class
            %
            % This function ultimately calls methodDocumentation.addMetadataFromDetailedDescription,
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
            arguments (Input)
                self ClassDocumentation
            end

            validMethods = [self.allMethodDocumentation(:).isHidden] == 0 & strcmp(string({self.allMethodDocumentation(:).access}),'public');
            if self.shouldExcludeAllSuperclasses == 1
                validMethods = validMethods & strcmp(string({self.allMethodDocumentation(:).declaringClassName}),self.name);
            else
                for iExcluded = 1:length(self.excludedSuperclasses)
                    validMethods = validMethods & ~strcmp(string({self.allMethodDocumentation(:).declaringClassName}),self.excludedSuperclasses{iExcluded});
                end
            end

            self.methodDocumentationNameMap = containers.Map;
            for index = find(validMethods)
                self.addMethodDocumentation(self.allMethodDocumentation(index));
            end

            % Use the detailed description to build the topic list
            mc = meta.class.fromName(self.name);
            self.rootTopic = Topic.topicsFromString(mc.DetailedDescription);
            ClassDocumentation.addPropertyAndMethodsToTopics(self.rootTopic,self.methodDocumentationNameMap);
        end

        function writeToFile(self)
            self.initializeMethodDocumentation();

            if ~exist(self.pathOfClassFolderAbsolute,'dir')
                mkdir(self.pathOfClassFolderAbsolute);
            end

            self.writeMarkdownForClassIndex();

            iPageNumber = 0;
            methodNames = keys(self.methodDocumentationNameMap);
            for i=1:length(methodNames)
                iPageNumber = iPageNumber+1;
                methodDocumentation = self.methodDocumentationNameMap(methodNames{i});
                methodDocumentation.writeToFile(iPageNumber)
            end
        end

        function writeMarkdownForClassIndex(self)
            arguments
                self
            end

            pathOfIndexFile = sprintf('%s/index.md',self.pathOfClassFolderAbsolute);

            fileID = fopen(pathOfIndexFile,'w');
            fprintf(fileID,'---\nlayout: default\ntitle: %s\nhas_children: false\nhas_toc: false\nmathjax: true\n',self.name);
            if ~isempty(self.parent)
                fprintf(fileID,'parent: %s\n',self.parent);
            end
            if ~isempty(self.grandparent)
                fprintf(fileID,'grand_parent: %s\n',self.grandparent);
            end
            if ~isempty(self.nav_order)
                fprintf(fileID,'nav_order: %d\n',self.nav_order);
            end
            fprintf(fileID,'---\n\n');

            fprintf(fileID,'#  %s\n',self.name);
            fprintf(fileID,'\n%s\n',self.shortDescription);

            fprintf(fileID,'\n\n---\n\n');

            if ~isempty(self.declaration)
                fprintf(fileID,'## Declaration\n');
                %fprintf(fileID,'```matlab\n%s\n```\n\n',declaration);
                % unfortunately we cannot directly put links into a markdown code block, so
                % we have to manually extract the url ourselves, and then write our own
                % html. so annoying!!
                mdurlExpression = '\[(?<name>[^\[]+)\]\((?<url>.*)\)';
                matchStr = regexpi(self.declaration,mdurlExpression,'names');
                if ~isempty(matchStr)
                    linkString = strcat('<a href="',matchStr.url,'" title="',matchStr.name,'">',matchStr.name,'</a>');
                    self.declaration = regexprep(self.declaration,mdurlExpression,linkString,'ignorecase');
                end
                prestring = '<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>';
                poststring = '</code></pre></div></div>';
                fprintf(fileID,'\n%s\n\n',strcat(prestring,strip(self.declaration),poststring));
            end

            if ~isempty(self.detailedDescription)
                fprintf(fileID,'## Overview\n%s\n',self.detailedDescription);
            end
            fprintf(fileID,'\n\n## Topics\n');

            % do not call the method directly, because we do not want to print the Root
            % category name
            for iSubtopic = 1:length(self.rootTopic.subtopics)
                ClassDocumentation.writeMarkdownForTopic(self.rootTopic.subtopics(iSubtopic),fileID,'',self.pathOfClassFolderRelative);
            end

            fprintf(fileID,'\n\n---');

            fclose(fileID);
        end


    end

    methods (Static)
        function declaration = declarationFromString(string)
            declarationExpression = '- declaration:(?<declaration>[^\r\n]+)(?:$|\n)';
            matchStr = regexpi(string,declarationExpression,'names');
            if ~isempty(matchStr)
                declaration = matchStr.declaration;
            else
                declaration = [];
            end
        end

        function trimmedString = trimDeclarationFromString(string)
            declarationExpression = '- declaration:(?<declaration>[^\r\n]+)(?:$|\n)';
            trimmedString  = regexprep(string,declarationExpression,'','ignorecase');
        end

        function classDocumentation = classDocumentationFromClassNames(nameArray,options)
            arguments
                nameArray 
                options.buildFolder % the folder where we are dumping everything on the local hard drive. This will become the *root* website folder
                options.websiteFolder % the folder relative to the root website folder
                options.parent = []
                options.grandparent = []
            end

            classDocumentation = ClassDocumentation.empty(length(nameArray),0);
            for iName=1:length(nameArray)
                classDocumentation(iName) = ClassDocumentation(nameArray{iName},buildFolder=options.buildFolder,websiteFolder=options.websiteFolder,parent=options.parent,grandparent=options.grandparent,nav_order=iName);
            end
        end

        function writeMarkdownForTopic(topic,fileID,indentLevel,websiteFolder)

            if isempty(topic.methodAnnotations) && isempty(topic.subtopics)
                return
            end
            fprintf(fileID,'%s+ %s\n',indentLevel,topic.name);
            for methodIndex = 1:length(topic.methodAnnotations)
                fprintf(fileID,'%s  + [`%s`](%s) ',indentLevel,topic.methodAnnotations(methodIndex).name,topic.methodAnnotations(methodIndex).pathOfFileOnWebsite);
                fprintf(fileID,'%s\n',topic.methodAnnotations(methodIndex).shortDescription);
            end
            for iSubtopic = 1:length(topic.subtopics)
                ClassDocumentation.writeMarkdownForTopic(topic.subtopics(iSubtopic),fileID,[indentLevel,'  '],websiteFolder);
            end

        end

        function addPropertyAndMethodsToTopics(rootTopic,metadataNameMap)
            % Adds methods and property metadata to the topic list
            %
            % - Parameter rootTopic: the rootTopic returned by topicsFromClass
            % - Parameter metadataNameMap: map with keys of method names and values of normalized metadata structs.
            arguments
                rootTopic Topic
                metadataNameMap containers.Map
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % Step 2) take the method/property/dimension/etc metadata sort those into
            % the appropriate topics.
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            methodNames = keys(metadataNameMap);
            otherTopic = Topic('Other');
            for i=1:length(methodNames)
                metadata = metadataNameMap(methodNames{i});

                % if the method has no topic, assign it to the 'other' topic we created
                if isempty(metadata.topic)
                    otherTopic.addMethod(metadata);
                    continue;
                end

                % if the topic does not exist, create it
                if isempty(rootTopic.subtopicWithName(metadata.topic))
                    rootTopic.addSubtopic(Topic(metadata.topic));
                end
                topic = rootTopic.subtopicWithName(metadata.topic);

                % if there is no subtopic, assign to the topic, and then exit
                if isempty(metadata.subtopic)
                    topic.addMethod(metadata);
                    continue;
                end

                % if the subtopic does not exist, create it
                if isempty(topic.subtopicWithName(metadata.subtopic))
                    topic.addSubtopic(Topic(metadata.subtopic));
                end
                subtopic = topic.subtopicWithName(metadata.subtopic);

                % if there is no subsubtopic, assign to the subtopic, and then exit
                if isempty(metadata.subsubtopic)
                    subtopic.addMethod(metadata);
                    continue;
                end

                % if the subsubtopic does not exist, create it
                if isempty(subtopic.subtopicWithName(metadata.subsubtopic))
                    subtopic.addSubtopic(Topic(metadata.subsubtopic));
                end
                subsubtopic = subtopic.subtopicWithName(metadata.subsubtopic);
                subsubtopic.addMethod(metadata);
            end

            % Make the 'Other' topic go at the end
            rootTopic.addSubtopic(otherTopic);
        end

        function methodDocumentation = methodDocumentationFromClass(className)
            methodDocumentation = MethodDocumentation.empty(0,0);
            classMetadata = meta.class.fromName(className);
            
            % add all the methods and properties
            for i=1:length(classMetadata.MethodList)
                mp = classMetadata.MethodList(i);
                metadata = MethodDocumentation(mp.Name);
                metadata.addMetadataFromMethodMetadata(mp);
                metadata.addMetadataFromDetailedDescription(mp.DetailedDescription);
                metadata.declaringClassName = classMetadata.Name;
                methodDocumentation(end+1) = metadata;
            end
            for i=1:length(classMetadata.PropertyList)
                mp = classMetadata.PropertyList(i);
                metadata = MethodDocumentation(mp.Name);
                metadata.addMetadataFromPropertyMetadata(mp);
                metadata.addMetadataFromDetailedDescription(mp.DetailedDescription);
                metadata.declaringClassName = classMetadata.Name;
                methodDocumentation(end+1) = metadata;
            end

            for iClass=1:length(classMetadata.SuperclassList)
                methodDocumentation = ClassDocumentation.annotateMethodDocumentationFromClassMetadata(methodDocumentation,classMetadata.SuperclassList(iClass));
            end
        end

        function methodDocumentation = annotateMethodDocumentationFromClassMetadata(methodDocumentation,classMetadata)
            subclassMethodNames = string({methodDocumentation(:).name});
            superclassMethodNames = string(cat(2,{classMetadata.MethodList(:).Name},{classMetadata.PropertyList(:).Name}));
            subclassIndices = ismember(subclassMethodNames,superclassMethodNames);
            for index = find(subclassIndices)
                methodDocumentation(index).declaringClassName = classMetadata.Name;
            end

            for iClass=1:length(classMetadata.SuperclassList)
                methodDocumentation = ClassDocumentation.annotateMethodDocumentationFromClassMetadata(methodDocumentation,classMetadata.SuperclassList(iClass));
            end
        end

    end
end