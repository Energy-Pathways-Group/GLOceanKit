classdef NetCDFComplexVariable < NetCDFVariable
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        realp NetCDFRealVariable
        imagp NetCDFRealVariable
    end

    properties (Dependent)
        value
    end

    methods
        function self = NetCDFComplexVariable(options)
            arguments
                options.real NetCDFRealVariable
                options.imag NetCDFRealVariable
                options.group (1,1) NetCDFGroup
                options.name {mustBeText}
                options.dimensions (:,1) NetCDFDimension
                options.attributes containers.Map = containers.Map(KeyType='char',ValueType='any');
                options.type string
            end
            if isfield(options,'real') && isfield(options,'imag')
                % this means we are initializing from file
                self.name = options.name;
                self.group = options.group;
                self.realp = options.real;
                self.imagp = options.imag;
                self.dimensions = self.realp.dimensions;
                self.type = self.realp.type;
                self.attributes = containers.Map(self.realp.attributes.keys,self.realp.attributes.values);
                self.attributes.remove(NetCDFVariable.GLNetCDFSchemaIsComplexKey);
                self.attributes.remove(NetCDFVariable.GLNetCDFSchemaIsRealPartKey);
                self.attributes.remove(NetCDFVariable.GLNetCDFSchemaIsImaginaryPartKey);
            else
                requiredFields = {'group','name','dimensions','type'};
                for iField=1:length(requiredFields)
                    if ~isfield(options,requiredFields{iField})
                        error('You must specify %s to initial a new variable.',requiredFields{iField});
                    end
                end

                self.group = options.group;
                self.dimensions = options.dimensions;
                self.name = options.name;
                self.type = options.type;
                self.attributes = options.attributes;

                realAttributes = containers.Map(KeyType='char',ValueType='any');
                imagAttributes = containers.Map(KeyType='char',ValueType='any');
                keyNames = self.attributes.keys;
                for iKey = 1:length(keyNames)
                    realAttributes(keyNames{iKey})=self.attributes(keyNames{iKey});
                    imagAttributes(keyNames{iKey})=self.attributes(keyNames{iKey});
                end

                realAttributes(NetCDFVariable.GLNetCDFSchemaIsComplexKey) = uint8(1);
                realAttributes(NetCDFVariable.GLNetCDFSchemaIsRealPartKey) = uint8(1);
                realAttributes(NetCDFVariable.GLNetCDFSchemaIsImaginaryPartKey) = uint8(0);
                self.realp =  NetCDFRealVariable(self.group,name=strcat(options.name,"_real"),dimensions=options.dimensions,attributes=realAttributes,type=options.type);

                
                imagAttributes(NetCDFVariable.GLNetCDFSchemaIsComplexKey) = uint8(1);
                imagAttributes(NetCDFVariable.GLNetCDFSchemaIsRealPartKey) = uint8(0);
                imagAttributes(NetCDFVariable.GLNetCDFSchemaIsImaginaryPartKey) = uint8(1);
                self.imagp =  NetCDFRealVariable(self.group,name=strcat(options.name,"_imag"),dimensions=options.dimensions,attributes=imagAttributes,type=options.type);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Add attributes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function addAttribute(self,key,value)
            self.attributes(key) = value;
            self.realp.addAttribute(key,value);
            self.imagp.addAttribute(key,value);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Read/write data
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function value = get.value(self)
            arguments
                self NetCDFComplexVariable
            end
            value = complex(self.realp.value,self.imagp.value);
        end

        function set.value(self,data)
            arguments
                self NetCDFComplexVariable {mustBeNonempty}
                data {mustBeNonempty}
            end
            self.realp.value = real(data);
            self.imagp.value = imag(data);
        end

        function value = valueAlongDimensionAtIndex(self,dimensionName,index)
            arguments
                self NetCDFComplexVariable {mustBeNonempty} 
                dimensionName char {mustBeNonempty}
                index {mustBeNonnegative}
            end
            value = complex(self.realp.valueAlongDimensionAtIndex(dimensionName,index),self.imagp.valueAlongDimensionAtIndex(dimensionName,index));
        end

        function setValueAlongDimensionAtIndex(self,data,dimensionName,index)
            % append new data to an existing variable
            %
            % concatenates data along a variable dimension (such as a time
            % dimension).
            %
            % ```matlab
            % variable.concatenateValueAlongDimensionAtIndex(data,dim,outputIndex);
            % ```
            %
            % - Topic: Working with variables
            % - Declaration: concatenateValueAlongDimensionAtIndex(data,dimension,index)
            % - Parameter data: variable data
            % - Parameter dimension: the variable dimension along which to concatenate
            % - Parameter index: index at which to write data
            arguments
                self NetCDFComplexVariable {mustBeNonempty} 
                data 
                dimensionName char {mustBeNonempty}
                index {mustBeNonnegative}
            end
            self.realp.setValueAlongDimensionAtIndex(real(data),dimensionName,index);
            self.imagp.setValueAlongDimensionAtIndex(imag(data),dimensionName,index);
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
                        if ~isempty(complexName)
                            imagName = strcat(complexName,"_imagp");
                        else
                            complexName = extractBefore(variable.name,"_real");
                            imagName = strcat(complexName,"_imag");
                        end                 
                        imagVariable = [];
                        for jVar=1:length(variables)
                            if strcmp(variables(jVar).name,imagName)
                                imagVariable = variables(jVar);
                            end
                        end
                        if ~isempty(imagVariable)
                            complexVariable = NetCDFComplexVariable(name=complexName,real=variable,imag=imagVariable,group=variable.group);
                            complexVariables(end+1) = complexVariable;
                        end
                    end
                end
            end
            complexVariables = reshape(complexVariables,[],1);
        end
    end
end