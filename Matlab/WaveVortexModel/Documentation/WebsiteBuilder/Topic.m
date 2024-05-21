classdef Topic < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        name char
        methodAnnotations
        subtopics
    end

    methods
        function self = Topic(name)
            self.name = name;
            self.subtopics = Topic.empty(0,0);
            self.methodAnnotations = MethodDocumentation.empty(0,0);
        end

        function addMethod(self,methodAnnotation)
            arguments
                self Topic
                methodAnnotation {mustBeNonempty}
            end
            self.methodAnnotations(end+1) = methodAnnotation;
            [~,indices] = sort([self.methodAnnotations.nav_order]);
            self.methodAnnotations = self.methodAnnotations(indices);
        end

        function addSubtopic(self,subtopic)
            arguments
                self Topic
                subtopic Topic
            end
            self.subtopics(end+1) = subtopic;
        end

        function subtopic = subtopicWithName(self,subtopicName)
            arguments (Input)
                self Topic
                subtopicName char
            end
            subtopic = [];
            for iSubtopic=1:length(self.subtopics)
                if strcmpi(self.subtopics(iSubtopic).name,subtopicName)
                    subtopic = self.subtopics(iSubtopic);
                end
            end
        end

        function description(self,options)
            arguments
                self Topic
                options.indent char = '';
            end
            fprintf('%s- %s\n',options.indent,self.name);
            for iMethod=1:length(self.methodAnnotations)
                fprintf('%s  - %s\n',options.indent,self.methodAnnotations(iMethod).name);
            end
            for iSubtopic=1:length(self.subtopics)
                self.subtopics(iSubtopic).description(indent=[options.indent,'  ']);
            end
        end
    end

    methods (Static)
        function detailedDescription = trimTopicsFromString(detailedDescription)
            % Remove any topic metadata from a string
            arguments
                detailedDescription
            end
            % extract topics and the detailed description (minus those topics)
            subsubtopicExpression = '- topic:([ \t]*)(?<topicName>[^—\r\n]+)—([ \t]*)(?<subtopicName>[^\r\n]+)—([ \t]*)(?<subsubtopicName>[^\r\n]+)(?:$|\n)';
            detailedDescription = regexprep(detailedDescription,subsubtopicExpression,'','ignorecase');
            subtopicExpression = '- topic:([ \t]*)(?<topicName>[^—\r\n]+)—([ \t]*)(?<subtopicName>[^\r\n]+)(?:$|\n)';
            detailedDescription = regexprep(detailedDescription,subtopicExpression,'','ignorecase');
            topicExpression = '- topic:([ \t]*)(?<topicName>[^\r\n]+)(?:$|\n)';
            detailedDescription = regexprep(detailedDescription,topicExpression,'','ignorecase');
        end

        function rootTopic = topicsFromString(detailedDescription)
            % Extracts topics and subtopic from a detailedDescription and creates a structure useful creating an ordered topic index
            %
            % - Parameter mc: the detailed description
            % - Returns detailedDescription: the detailed description without the topic/subtopic metadata
            % - Returns rootTopic: a Topic instance containing the topics and subtopics defined by the class
            arguments
                detailedDescription
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % Step 1) extract a list of topics and subtopics from the class's
            % detailedDescription, and then sort those into useful structures.
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % extract topics and the detailed description (minus those topics)
            subsubtopicExpression = '- topic:([ \t]*)(?<topicName>[^—\r\n]+)—([ \t]*)(?<subtopicName>[^\r\n]+)—([ \t]*)(?<subsubtopicName>[^\r\n]+)(?:$|\n)';
            classDefinedSubsubTopics = regexpi(detailedDescription,subsubtopicExpression,'names');
            detailedDescription = regexprep(detailedDescription,subsubtopicExpression,'','ignorecase');
            subtopicExpression = '- topic:([ \t]*)(?<topicName>[^—\r\n]+)—([ \t]*)(?<subtopicName>[^\r\n]+)(?:$|\n)';
            classDefinedSubTopics = regexpi(detailedDescription,subtopicExpression,'names');
            detailedDescription = regexprep(detailedDescription,subtopicExpression,'','ignorecase');
            topicExpression = '- topic:([ \t]*)(?<topicName>[^\r\n]+)(?:$|\n)';
            classDefinedTopics = regexpi(detailedDescription,topicExpression,'names');

            rootTopic = Topic('Root');
            for iTopic=1:length(classDefinedTopics)
                rootTopic.addSubtopic(Topic(strtrim(classDefinedTopics(iTopic).topicName)));
            end

            for iSubtopic=1:length(classDefinedSubTopics)
                topic = rootTopic.subtopicWithName(strtrim(classDefinedSubTopics(iSubtopic).topicName));
                topic.addSubtopic(Topic(strtrim(classDefinedSubTopics(iSubtopic).subtopicName)));
            end

            for iSubtopic=1:length(classDefinedSubsubTopics)
                topic = rootTopic.subtopicWithName(strtrim(classDefinedSubsubTopics(iSubtopic).topicName));
                subtopic = topic.subtopicWithName(strtrim(classDefinedSubsubTopics(iSubtopic).subtopicName));
                subtopic.addSubtopic(Topic(strtrim(classDefinedSubsubTopics(iSubtopic).subsubtopicName)));
            end
        end
    end
end