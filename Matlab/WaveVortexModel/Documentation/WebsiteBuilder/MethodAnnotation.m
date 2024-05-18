classdef MethodAnnotation < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        name string
        className
        parameters
        returns
        detailedDescription =[]
        shortDescription = []
        declaration = []
        topic = []
        subtopic = []
        subsubtopic = []
        nav_order = Inf
        functionType

        dimensions
        units
        isComplex
    end

    methods
        function self = MethodAnnotation(name)
            self.name = name;
        end

        function self = addMetadataFromDetailedDescription(self,detailedDescription)
            if isempty(detailedDescription)
                return;
            end

            % Check out https://regexr.com for testing these regex.
            topicExpression = '- topic:([ \t]*)(?<topic>[^\r\n]+)(?:$|\n)';
            subtopicExpression = '- topic:([ \t]*)(?<topic>[^—\r\n]+)—([ \t]*)(?<subtopic>[^\r\n]+)(?:$|\n)';
            subsubtopicExpression = '- topic:([ \t]*)(?<topic>[^—\r\n]+)—([ \t]*)(?<subtopic>[^\r\n]+)—([ \t]*)(?<subsubtopic>[^\r\n]+)(?:$|\n)';
            declarationExpression = '- declaration:(?<declaration>[^\r\n]+)(?:$|\n)';
            parameterExpression = '- parameter (?<name>[^:]+):(?<description>[^\r\n]+)(?:$|\n)';
            returnsExpression = '- returns (?<name>[^:]+):(?<description>[^\r\n]+)(?:$|\n)';
            navOrderExpression = '- nav_order:([ \t]*)(?<nav_order>[^\r\n]+)(?:$|\n)';
            leadingWhitespaceExpression = '^[ \t]+';

            % Capture the subsubtopic annotation, then remove it
            matchStr = regexpi(detailedDescription,subsubtopicExpression,'names');
            detailedDescription = regexprep(detailedDescription,subsubtopicExpression,'','ignorecase');
            if ~isempty(matchStr)
                self.subsubtopic = strtrim(matchStr.subsubtopic);
                self.subtopic = strtrim(matchStr.subtopic);
                self.topic = strtrim(matchStr.topic);
            end

            % Capture the subtopic annotation, then remove it
            matchStr = regexpi(detailedDescription,subtopicExpression,'names');
            detailedDescription = regexprep(detailedDescription,subtopicExpression,'','ignorecase');
            if ~isempty(matchStr)
                self.subtopic = strtrim(matchStr.subtopic);
                self.topic = strtrim(matchStr.topic);
            end

            % Capture the topic annotation, then remove it
            matchStr = regexpi(detailedDescription,topicExpression,'names');
            detailedDescription = regexprep(detailedDescription,topicExpression,'','ignorecase');
            if ~isempty(matchStr)
                self.topic = strtrim(matchStr.topic);
            end

            % Capture all parameters, then remove the annotations
            self.parameters = regexpi(detailedDescription,parameterExpression,'names');
            detailedDescription = regexprep(detailedDescription,parameterExpression,'','ignorecase');

            % Capture all returns, then remove the annotations
            self.returns = regexpi(detailedDescription,returnsExpression,'names');
            detailedDescription = regexprep(detailedDescription,returnsExpression,'','ignorecase');


            % Capture any declarations made, then remove the annotation
            matchStr = regexpi(detailedDescription,declarationExpression,'names');
            detailedDescription = regexprep(detailedDescription,declarationExpression,'','ignorecase');
            if ~isempty(matchStr)
                self.declaration = matchStr.declaration;
            end

            matchStr = regexpi(detailedDescription,navOrderExpression,'names');
            detailedDescription = regexprep(detailedDescription,navOrderExpression,'','ignorecase');
            if ~isempty(matchStr)
                self.nav_order = str2double(matchStr.nav_order);
            end


            self.detailedDescription = regexprep(detailedDescription,leadingWhitespaceExpression,'');
        end
    end
end