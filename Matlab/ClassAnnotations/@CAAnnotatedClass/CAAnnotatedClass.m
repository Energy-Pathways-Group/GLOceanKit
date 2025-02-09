classdef CAAnnotatedClass < handle
    %CAAnnotatedClass A class with metadata describing its properties
    %
    % An `annotated class' has metadata describing its properties as either
    % dimensions or properties (which may have dimensions). This metadata
    % is used to read/write NetCDF files, and to create website
    % documentation.
    %
    % Annotations should specify the netcdf group they belong to.
    % Annotations should specify whether the type is a function_handle
    %
    % A subclass must do two things:
    % 1) implement the required methods (perhaps by simply calling the
    % static methods) and,
    % 2) initialize the superclass during init, self@CAAnnotatedClass()
    % 
    properties (Access=private)
        propertyAnnotationNameMap
        dimensionAnnotationNameMap
    end

    events
        propertyAnnotationsDidChange
    end

    methods
        function self = CAAnnotatedClass()
            % Matlab dictionaries do not allow a subclass type; the
            % work-around is to use a cell.
            self.dimensionAnnotationNameMap = configureDictionary("string","cell");
            self.propertyAnnotationNameMap = configureDictionary("string","cell");

            self.addPropertyAnnotation(feval(strcat(class(self),'.classDefinedPropertyAnnotations')));
        end
        
        vars = requiredProperties(self);

        names = annotatedDimensionNames(self)
        val = dimensionAnnotationWithName(self,name)
        function dimAnnotations = dimensionAnnotations(self)
            arguments (Input)
                self CAAnnotatedClass
            end
            arguments (Output)
                dimAnnotations CADimensionProperty
            end
            dimAnnotations = [self.dimensionAnnotationNameMap{self.dimensionAnnotationNameMap.keys}];
        end

        addPropertyAnnotation(self,variableAnnotation)
        removePropertyAnnotation(self,variableAnnotation)
        val = propertyAnnotationWithName(self,name)
        names = annotatedPropertyNames(self)
        function propAnnotations = propertyAnnotations(self)
            arguments (Input)
                self CAAnnotatedClass
            end
            arguments (Output)
                propAnnotations CAPropertyAnnotation
            end
            propAnnotations = [self.propertyAnnotationNameMap{self.propertyAnnotationNameMap.keys}];
        end

        ncfile = writeToFile(self,path,variables,options)
        
        function flag = isequal(self,other)
            arguments
                self CAAnnotatedClass
                other CAAnnotatedClass
            end
            flag = isequal(class(self),class(other));
            names = self.requiredProperties();
            for name = names
                if isa(self.(name{1}),"function_handle")
                    flag = flag & isequal(func2str(self.(name{1})), func2str(other.(name{1})));
                elseif isa(self.(name{1}),"CAAnnotatedClass")
                    myObjs = self.(name{1});
                    otherObjs = other.(name{1});
                    if length(myObjs) == length(otherObjs)
                        for iObj=1:length(myObjs)
                            flag = flag & isequal(myObjs(iObj),otherObjs(iObj));
                        end
                    else
                        flag = false;
                    end
                else
                    flag = flag & isequal(self.(name{1}), other.(name{1}));
                end
            end
        end
    end

    methods (Static)
        propertyAnnotations = classDefinedPropertyAnnotations()
        atc = annotatedClassFromFile(path)

        function atc = annotatedClassFromGroup(group)
            className = group.attributes('AnnotatedClass');
            requiredProperties = feval(strcat(className,'.classRequiredPropertyNames'));
            var = CAAnnotatedClass.propertyValuesFromGroup(group,requiredProperties);
            varCell = namedargs2cell(var);
            atc = feval(className,varCell{:});
        end

        function var = requiredPropertiesFromGroup(group,options)
            arguments (Input)
                group NetCDFGroup
                options.className string
            end
            arguments (Output)
                var struct
            end
            if isfield(options,'className')
                className = options.className;
            else
                if isKey(group.attributes,'AnnotatedClass')
                    className = group.attributes('AnnotatedClass');
                else
                    error('Unable to find the attribute AnnotatedClass');
                end
            end
            requiredProperties = feval(strcat(className,'.classRequiredPropertyNames'));
            var = CAAnnotatedClass.propertyValuesFromGroup(group,requiredProperties);
        end

        function var = propertyValuesFromGroup(group,propertyNames,options)
            arguments (Input)
                group NetCDFGroup
                propertyNames {mustBeA(propertyNames,"cell")} = {}
                options.shouldIgnoreMissingProperties logical = false
            end
            arguments (Output)
                var struct
            end
            for iVar = 1:length(propertyNames)
                name = propertyNames{iVar};
                if group.hasVariableWithName(name) == true
                    var.(name) = group.readVariables(name);
                    continue;
                end
                
                targetGroupName = join( [string(group.name),string(name)],"-");
                if group.parentGroup.hasGroupWithName(targetGroupName) == true
                    targetGroup = group.parentGroup.groupWithName(targetGroupName);
                    if isempty(targetGroup.groups)
                        if isKey(targetGroup.attributes,'AnnotatedClass')
                            var.(name) = CAAnnotatedClass.annotatedClassFromGroup(targetGroup);
                            continue;
                        end
                    else
                        for iObj=1:length(targetGroup.groups)
                            subGroupName = join( [string(name),string(iObj)],"-");
                            if targetGroup.hasGroupWithName(subGroupName)
                                subGroup = targetGroup.groupWithName(subGroupName);
                                if isKey(subGroup.attributes,'AnnotatedClass')
                                    tmp(iObj) = CAAnnotatedClass.annotatedClassFromGroup(subGroup);
                                end
                            end
                        end
                        var.(name) = tmp;
                        continue;
                    end
                end

                if options.shouldIgnoreMissingProperties == true
                    continue;
                else
                    error('Unable to find the property %s',name);
                end
            end
        end

        function [canInit, errorString] = canInitializeDirectlyFromGroup(group,requiredProperties)
            canInit = false;
            errorString = "";
            hasVar = group.hasVariableWithName(requiredProperties{:});
            if ~all(hasVar)
                errorString = "The NetCDF file is missing required variables: " + join(string(requiredProperties(~hasVar)),', ');
                return;
            end

            canInit = true;
        end

        function throwErrorIfMissingProperties(group,requiredProperties)
            arguments (Input)
                group NetCDFGroup
                requiredProperties {mustBeA(requiredProperties,"cell")} = {}
            end
            [canInit, errorString] = CAAnnotatedClass.canInitializeDirectlyFromGroup(group,requiredProperties);
            if ~canInit
                error(errorString);
            end
        end

        function ncfile = writeToPath(ac,path,options)
            arguments (Input)
                ac CAAnnotatedClass {mustBeNonempty}
                path char {mustBeNonempty}
                options.shouldOverwriteExisting double {mustBeMember(options.shouldOverwriteExisting,[0 1])} = 0
                options.propertyAnnotations CAPropertyAnnotation = CAPropertyAnnotation.empty(0,0)
                options.attributes = configureDictionary("string","string")
            end

            arguments (Output)
                ncfile NetCDFFile
            end

            if options.shouldOverwriteExisting == 1
                if isfile(path)
                    delete(path);
                end
            else
                if isfile(path)
                    error('A file already exists with that name.')
                end
            end
            ncfile = NetCDFFile(path);
            CAAnnotatedClass.writeToGroup(ac,ncfile,options.propertyAnnotations,options.attributes);
        end

        function writeToGroup(ac,targetGroup,propertyAnnotations,attributes,options)
            arguments
                ac CAAnnotatedClass
                targetGroup NetCDFGroup
                propertyAnnotations CAPropertyAnnotation = CAPropertyAnnotation.empty(0,0)
                attributes = configureDictionary("string","string")
                options.shouldAlwaysUseSubgroup logical = false
            end
            targetGroup.addAttribute('AnnotatedClass',class(ac));

            % if one or more of the properties are other CAAnnotatedClass
            % objects, then we have to stuff this group into its own
            % subgroup to prevent dimensions from cascading down.
            requiresSubgroup = false | options.shouldAlwaysUseSubgroup;
            for i=1:length(propertyAnnotations)
                if isa(propertyAnnotations(i),'CAObjectProperty')
                    requiresSubgroup = true;
                end
            end
            if requiresSubgroup
                group = targetGroup.addGroup(class(ac));
                group.addAttribute('AnnotatedClass',class(ac));
            else
                group = targetGroup;
            end

            attributeNames = keys(attributes);
            for iKey=1:length(attributeNames)
                group.addAttribute(attributeNames{iKey},attributes(attributeNames{iKey}));
            end

            for i=1:length(propertyAnnotations)
                if isa(propertyAnnotations(i),'CADimensionProperty')
                    dimAnnotation = propertyAnnotations(i);
                    dimAnnotation.attributes('units') = dimAnnotation.units;
                    dimAnnotation.attributes('long_name') = dimAnnotation.description;
                    group.addDimension(dimAnnotation.name,ac.(dimAnnotation.name),attributes=dimAnnotation.attributes);
                end
            end

            for i=1:length(propertyAnnotations)
                if isa(propertyAnnotations(i),'CAFunctionProperty')
                    group.addFunctionHandle(propertyAnnotations(i).name,ac.(propertyAnnotations(i).name),attributes=propertyAnnotations(i).attributes);
                elseif isa(propertyAnnotations(i),'CAObjectProperty')
                    obj = ac.(propertyAnnotations(i).name);
                    assert( isa(obj,'CAAnnotatedClass'),'The object property %s is not a subclass of CAAnnotatedClass, and thus cannot be added to file',propertyAnnotations(i).name);
                    objGroup = targetGroup.addGroup(join( [string(class(ac)),string(propertyAnnotations(i).name)],"-"));
                    if isscalar(obj)
                        objPropertyAnnotations = obj.propertyAnnotationWithName(obj.requiredProperties);
                        CAAnnotatedClass.writeToGroup(obj,objGroup,objPropertyAnnotations);
                    else
                        for iObj=1:length(obj)
                            sGroup = objGroup.addGroup(join( [string(propertyAnnotations(i).name),string(iObj)],"-"));
                            objPropertyAnnotations = obj(iObj).propertyAnnotationWithName(obj(iObj).requiredProperties);
                            CAAnnotatedClass.writeToGroup(obj(iObj),sGroup,objPropertyAnnotations);
                        end
                    end
                elseif isa(propertyAnnotations(i),'CANumericProperty')
                    propAttributes = propertyAnnotations(i).attributes;
                    if ~isempty(propertyAnnotations(i).units)
                        propAttributes('units') = propertyAnnotations(i).units;
                    end
                    if ~isempty(propertyAnnotations(i).description)
                        propAttributes('long_name') = propertyAnnotations(i).description;
                    end
                    group.addVariable(propertyAnnotations(i).name,propertyAnnotations(i).dimensions,ac.(propertyAnnotations(i).name),isComplex=propertyAnnotations(i).isComplex,attributes=propAttributes);
                end
            end

            % dimensionAnnotationDictionary = dictionary(string({options.dimAnnotations.name}),num2cell(options.dimAnnotations));
        end
    end
end