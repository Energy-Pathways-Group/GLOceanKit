classdef NetCDFComplexVariable < NetCDFVariable
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        real
        imag
    end

    properties (Dependent)
        value
    end

    methods
        function self = NetCDFComplexVariable(options)
            arguments
                options.real NetCDFVariable = []
                options.imag NetCDFVariable = []
                options.group (1,1) NetCDFGroup
                options.name string
                options.dimensions (:,1) NetCDFDimension
                options.attributes containers.Map = containers.Map(KeyType='char',ValueType='any');
                options.type string
            end
            if ~isempty(options.real) && ~isempty(options.imag)
                % this means we are initializing from file
                self.name = options.name;
                self.real = options.real;
                self.imag = options.imag;
                self.dimensions = self.real.dimensions;
                self.type = self.real.type;
                self.attributes = containers.Map(self.real.attributes.keys,self.real.attributes.values);
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

                % Now initialize!!, e.g.,
                realAttributes = containers.Map(self.attributes.keys,self.attributes.values);
                realAttributes(NetCDFVariable.GLNetCDFSchemaIsComplexKey) = uint8(1);
                realAttributes(NetCDFVariable.GLNetCDFSchemaIsRealPartKey) = uint8(1);
                realAttributes(NetCDFVariable.GLNetCDFSchemaIsImaginaryPartKey) = uint8(0);
                self.real =  NetCDFVariable(self,name=strcat(name,"_real"),dimensions=dims,attributes=realAttributes,type=self.type);

                imagAttributes = containers.Map(self.attributes.keys,self.attributes.values);
                imagAttributes(NetCDFVariable.GLNetCDFSchemaIsComplexKey) = uint8(1);
                imagAttributes(NetCDFVariable.GLNetCDFSchemaIsRealPartKey) = uint8(0);
                imagAttributes(NetCDFVariable.GLNetCDFSchemaIsImaginaryPartKey) = uint8(1);
                self.imag =  NetCDFVariable(self,name=strcat(name,"_imag"),dimensions=dims,attributes=imagAttributes,type=ncType);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Add attributes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function addAttribute(self,key,value)
            self.attributes(key) = value;
            self.real.addAttribute(key,value);
            self.imag.addAttribute(key,value);
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
            value = complex(self.real.value,self.imag.value);
        end

        function set.value(self,data)
            arguments
                self NetCDFComplexVariable {mustBeNonempty}
                data {mustBeNonempty}
            end
            self.real.value = real(data);
            self.imag.value = imag(data);
        end

        function value = valueAlongDimensionAtIndex(self,dimensionName,index)
            arguments
                self NetCDFComplexVariable {mustBeNonempty} 
                dimensionName char {mustBeNonempty}
                index {mustBeNonnegative}
            end
            value = complex(self.real.valueAlongDimensionAtIndex(dimensionName,index),self.imag.valueAlongDimensionAtIndex(dimensionName,index));
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
            self.real.setValueAlongDimensionAtIndex(real(data),dimensionName,index);
            self.imag.setValueAlongDimensionAtIndex(imag(data),dimensionName,index);
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
                            complexVariable = NetCDFComplexVariable(name=complexName,real=variable,imag=imagVariable);
                            complexVariables(end+1) = complexVariable;
                        end
                    end
                end
            end
            complexVariables = reshape(complexVariables,[],1);
        end
    end
end