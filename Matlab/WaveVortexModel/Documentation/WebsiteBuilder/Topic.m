classdef Topic < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        name char
        methodNames
        methodOrder
        subtopics
    end

    methods
        function self = Topic(name)
            self.name = name;
            self.subtopics = Topic.empty(0,0);
            self.methodNames = cell(0);
            self.methodOrder = cell(0);
        end

        function addMethod(self,metadata)
            arguments
                self Topic
                metadata {mustBeNonempty}
            end
            self.methodNames{end+1} = metadata.name;
            if isfield(metadata,'nav_order')
                self.methodOrder{end+1} = metadata.nav_order;
            else
                self.methodOrder{end+1} = Inf;
            end
            [self.methodOrder,indices] = sortrows(self.methodOrder);
            self.methodNames = self.methodNames(indices,:);
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
            for iMethod=1:length(self.methodNames)
                fprintf('%s  - %s\n',options.indent,self.methodNames{iMethod});
            end
            for iSubtopic=1:length(self.subtopics)
                self.subtopics(iSubtopic).description(indent=[options.indent,'  ']);
            end
        end
    end
end