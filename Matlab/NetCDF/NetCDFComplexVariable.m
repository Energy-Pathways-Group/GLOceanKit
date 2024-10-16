classdef NetCDFComplexVariable < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        name
        realVar
        imagVar
    end

    methods
        function self = NetCDFComplexVariable(name,realVar,imagVar)
            self.name = name;
            self.realVar = realVar;
            self.imagVar = imagVar;
        end

    end

    methods (Static)
        function complexVariables = complexVariablesFromVariables(variables)
            complexVariables = NetCDFComplexVariable.empty(0,0);
            for iVar=1:length(variables)
                variable = variables(iVar);
                if isKey(variable.attributes,{variable.GLNetCDFSchemaIsComplexKey,variable.GLNetCDFSchemaIsRealPartKey})
                    if variable.attributes(variable.GLNetCDFSchemaIsComplexKey) == 1 && variable.attributes(variable.GLNetCDFSchemaIsRealPartKey) == 1
                        complexName = extractBefore(variable.name,"_realp");
                        imagName = strcat(complexName,"_imagp");
                        imagVariable = [];
                        for jVar=1:length(variables)
                            if strcmp(variables(jVar).name,imagName)
                                imagVariable = variables(jVar);
                            end
                        end
                        if ~isempty(imagVariable)
                            complexVariable = NetCDFComplexVariable(complexName,variable,imagVariable);
                            complexVariables(end+1) = complexVariable;
                        end
                    end
                end
            end
            complexVariables = reshape(complexVariables,[],1);
        end
    end
end