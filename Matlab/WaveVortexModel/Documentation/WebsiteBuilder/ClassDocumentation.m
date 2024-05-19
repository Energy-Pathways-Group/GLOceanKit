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

        methodDocumentationNameMap

        shortDescription
        detailedDescription
        declaration = []
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
            self.detailedDescription = mc.DetailedDescription;
            self.shortDescription = mc.Description;

            declarationExpression = '- declaration:(?<declaration>[^\r\n]+)(?:$|\n)';
            matchStr = regexpi(self.detailedDescription,declarationExpression,'names');
            if ~isempty(matchStr)
                self.declaration = matchStr.declaration;
            end
            self.detailedDescription  = regexprep(self.detailedDescription,declarationExpression,'','ignorecase');
     
            % Use the detailed description to build the topic list
            [self.rootTopic,self.detailedDescription] = Topic.topicsFromDetailedDescription(self.detailedDescription);

            self.methodDocumentationNameMap = ClassDocumentation.metadataFromClassPropertiesAndMethods(mc);
            ClassDocumentation.addPropertyAndMethodsToTopics(self.rootTopic,self.methodDocumentationNameMap);
        end

        function writeToFile(self)
            self.writeMarkdownForClass();

            iPageNumber = 0;
            methodNames = keys(self.methodDocumentationNameMap);
            for i=1:length(methodNames)
                iPageNumber = iPageNumber+1;
                methodDocumentation = self.methodDocumentationNameMap(methodNames{i});
                methodDocumentation.writeToFile(sprintf('%s/%s.md',self.pathOfClassFolderAbsolute,lower(methodDocumentation.name)),iPageNumber)
            end
        end

        function writeMarkdownForClass(self)
            arguments
                self
            end

            if ~exist(self.pathOfClassFolderAbsolute,'dir')
                mkdir(self.pathOfClassFolderAbsolute);
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
                fprintf(fileID,'%s  + [`%s`](/%s/%s.html) ',indentLevel,topic.methodAnnotations(methodIndex).name,websiteFolder,lower(topic.methodAnnotations(methodIndex).name));
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
            arguments (Input)
                mc meta.class
            end
            arguments (Output)
                metadataNameMap containers.Map
            end

            metadataNameMap = containers.Map;
            for i=1:length(mc.MethodList)
                metadata = ClassDocumentation.extractMethodMetadata(mc.MethodList(i),mc.Name);
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
                metadata = ClassDocumentation.extractMethodMetadata(mc.PropertyList(i),mc.Name);
                if ~isempty(metadata)
                    metadata.functionType = FunctionType.instanceProperty;
                    metadataNameMap(metadata.name) = metadata;
                end
            end
        end

        function metadata = extractMethodMetadata(mp,className)
            % Extract documentation from method or property (mp) metadata.
            metadata = [];

            % Don't create documentation if this is a method defined in the superclass
            % This initial check does not work, because if the subclass re-defines a
            % method, then it counts the subclass as the "DefiningClass". But, for
            % documentation purposes, we really don't want that.
            if ~strcmp(mp.DefiningClass.Name,className)
                return;
            end
            if ClassDocumentation.isMethodDefinedInSuperclass(mp)
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

            metadata = MethodDocumentation(mp.Name);
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

    end
end