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
        pathOfClassFolderOnWebsite
        pathOfClassFolderOnHardDrive
        parent
        grandparent
        nav_order
        rootTopic

        allMethodDocumentation
        shouldExcludeAllSuperclasses

        shortDescription
        detailedDescription
        declaration = []
    end

    properties (SetObservable, AbortSet, Access = public)
        excludedSuperclasses
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
                options.excludedSuperclasses = {'handle'};
            end
            self.name = name;
            self.parent = options.parent;
            self.grandparent = options.grandparent;
            self.nav_order = options.nav_order;
            self.excludedSuperclasses = options.excludedSuperclasses;
            if isequal(options.excludedSuperclasses,{'handle'})
                self.shouldExcludeAllSuperclasses = 1;
            else
                self.shouldExcludeAllSuperclasses = 0;
            end

            % relative path to a folder for all the class contents
            self.pathOfClassFolderOnWebsite = fullfile(options.websiteFolder,lower(self.name));

            % Make a folder for all the class contents
            self.pathOfClassFolderOnHardDrive = fullfile(options.buildFolder,self.pathOfClassFolderOnWebsite);

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
            % Add new method documentation beyond what was automatically
            % extracted from the class metadata---overwriting any existing
            % documentation with the same (case-sensitive) name.
            existsAtIndex = strcmp(string({self.allMethodDocumentation(:).name}),methodDocumentation.name);
            if any(existsAtIndex)
                index = find(existsAtIndex);
                self.allMethodDocumentation(index) = methodDocumentation;
            else
                self.allMethodDocumentation(end+1) = methodDocumentation;
            end
        end

        function writeToFile(self)
            % First identify which methods we want to include
            validMethods = [self.allMethodDocumentation(:).isHidden] == 0 & strcmp(string({self.allMethodDocumentation(:).access}),'public');
            if self.shouldExcludeAllSuperclasses == 1
                validMethods = validMethods & arrayfun(@(a) a.isDeclaredInClass(self.name),self.allMethodDocumentation);
            else
                for iExcluded = 1:length(self.excludedSuperclasses)
                    validMethods = validMethods & ~arrayfun(@(a) a.isDeclaredInClass(self.excludedSuperclasses{iExcluded}),self.allMethodDocumentation);
                end
            end

            % Extract those methods and assign a path to write.
            methodDocumentationNameMap = containers.Map;
            for index = find(validMethods)
                methodDocumentation = self.allMethodDocumentation(index);
                if ismember(lower(methodDocumentation.name),lower(keys(methodDocumentationNameMap))) && ~isKey(methodDocumentationNameMap,methodDocumentation.name)
                    methodDocumentation.pathOfOutputFile = fullfile(self.pathOfClassFolderOnHardDrive,sprintf('%s_.md',lower(methodDocumentation.name)));
                    methodDocumentation.pathOfFileOnWebsite = fullfile("/",self.pathOfClassFolderOnWebsite,sprintf('%s_.html',lower(methodDocumentation.name)));
                else
                    methodDocumentation.pathOfOutputFile = fullfile(self.pathOfClassFolderOnHardDrive,sprintf('%s.md',lower(methodDocumentation.name)));
                    methodDocumentation.pathOfFileOnWebsite = fullfile("/",self.pathOfClassFolderOnWebsite,sprintf('%s.html',lower(methodDocumentation.name)));
                end
                methodDocumentationNameMap(methodDocumentation.name) = methodDocumentation;
            end

            % Use the *class* detailed description to build the topic list
            mc = meta.class.fromName(self.name);
            self.rootTopic = Topic.topicsFromString(mc.DetailedDescription);
            % Then use the *method* information to add new topics and
            % assign methods to existing topics.
            ClassDocumentation.addPropertyAndMethodsToTopics(self.rootTopic,methodDocumentationNameMap);

            if ~exist(self.pathOfClassFolderOnHardDrive,'dir')
                mkdir(self.pathOfClassFolderOnHardDrive);
            end

            self.writeMarkdownForClassIndex();

            iPageNumber = 0;
            methodNames = keys(methodDocumentationNameMap);
            for i=1:length(methodNames)
                iPageNumber = iPageNumber+1;
                methodDocumentation = methodDocumentationNameMap(methodNames{i});
                methodDocumentation.writeToFile(self.name,iPageNumber)
            end
        end

        function writeMarkdownForClassIndex(self)
            arguments
                self
            end

            pathOfIndexFile = sprintf('%s/index.md',self.pathOfClassFolderOnHardDrive);

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
                ClassDocumentation.writeMarkdownForTopic(self.rootTopic.subtopics(iSubtopic),fileID,'',self.pathOfClassFolderOnWebsite);
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

            classDocumentation = ClassDocumentation.empty(length(nameArray),0);
            for iName=1:length(nameArray)
                classDocumentation(iName) = ClassDocumentation(nameArray{iName},buildFolder=options.buildFolder,websiteFolder=options.websiteFolder,parent=options.parent,grandparent=options.grandparent,nav_order=iName,excludedSuperclasses = options.excludedSuperclasses);
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
                if isequal(mp.DefiningClass.Name,classMetadata.Name)
                    metadata.addDeclaringClass(classMetadata.Name);
                end
                methodDocumentation(end+1) = metadata;
            end
            for i=1:length(classMetadata.PropertyList)
                mp = classMetadata.PropertyList(i);
                metadata = MethodDocumentation(mp.Name);
                metadata.addMetadataFromPropertyMetadata(mp);
                metadata.addMetadataFromDetailedDescription(mp.DetailedDescription);
                if isequal(mp.DefiningClass.Name,classMetadata.Name)
                    metadata.addDeclaringClass(classMetadata.Name);
                end
                methodDocumentation(end+1) = metadata;
            end

            for iClass=1:length(classMetadata.SuperclassList)
                methodDocumentation = ClassDocumentation.annotateMethodDocumentationFromClassMetadata(methodDocumentation,classMetadata.SuperclassList(iClass));
            end
        end

        function methodDocumentation = annotateMethodDocumentationFromClassMetadata(methodDocumentation,classMetadata)
            subclassMethodNames = string({methodDocumentation(:).name});

            % Which superclass methods have a defining class that is this
            % class
            superclassMethodNames = string(cat(2,{classMetadata.MethodList(:).Name},{classMetadata.PropertyList(:).Name}));
            a=cat(2,[classMetadata.MethodList(:).DefiningClass],[classMetadata.PropertyList(:).DefiningClass]);
            definingClassNames = string({a.Name});
            superclassMethodNames=superclassMethodNames(ismember(definingClassNames,classMetadata.Name));
            subclassIndices = ismember(subclassMethodNames,superclassMethodNames);
            for index = find(subclassIndices)
                methodDocumentation(index).addDeclaringClass(classMetadata.Name);
            end

            for iClass=1:length(classMetadata.SuperclassList)
                methodDocumentation = ClassDocumentation.annotateMethodDocumentationFromClassMetadata(methodDocumentation,classMetadata.SuperclassList(iClass));
            end
        end

    end
end