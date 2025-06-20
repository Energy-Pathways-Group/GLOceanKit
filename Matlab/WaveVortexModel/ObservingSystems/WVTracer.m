classdef WVTracer < WVObservingSystem
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties (GetAccess=public, SetAccess=protected)
        isXYOnly
        phi
        absTolerance
        shouldAntialias
    end

    methods
        function self = WVTracer(model,options)
            %create a new observing system
            %
            % This class is intended to be subclassed, so it generally
            % assumed that this initialization will not be called directly.
            %
            % - Topic: Initialization
            % - Declaration: self = WVObservingSystem(model,name)
            % - Parameter model: the WVModel instance
            % - Parameter name: name of the observing system
            % - Returns self: a new instance of WVObservingSystem
            arguments
                model WVModel
                options.name {mustBeText}
                options.phi double
                options.isXYOnly (1,1) logical = false
                options.absTolerance = 1e-5; % set assuming phi is normalized to a max value of 1
                options.shouldAntialias = true;
            end
            self@WVObservingSystem(model,options.name);
            if ~isfield(options,'phi')
                error('You must specify the initial tracer field phi');
            end

            self.isXYOnly = options.isXYOnly;
            self.phi = options.phi;
            self.absTolerance = options.absTolerance;
            self.shouldAntialias = options.shouldAntialias;

            self.nFluxComponents = 1;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Integrated variables
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function nArray = lengthOfFluxComponents(self)
            % return an array containing the numel of each flux component.
            nArray = numel(self.phi);
        end

        function Y0 = absErrorTolerance(self)
            Y0 = {self.absTolerance};
        end

        function Y0 = initialConditions(self)
            Y0 = {self.phi};
        end

        function F = fluxAtTime(self,t,y0)
            self.updateIntegratorValues(t,y0);
            if self.isXYOnly
                F_phi = -self.wvt.u .* self.wvt.diffX(self.phi) - self.wvt.v .* self.wvt.diffY(self.phi);
            else
                F_phi = -self.wvt.u .* self.wvt.diffX(self.phi) - self.wvt.v .* self.wvt.diffY(self.phi) - self.wvt.w .* self.wvt.diffZF(self.phi);
            end
            if self.shouldAntialias
                F_phi = self.wvt.transformToSpatialDomainWithFourier(self.wvt.transformFromSpatialDomainWithFourier(F_phi));
            end
            F = {F_phi};
        end

        function updateIntegratorValues(self,t,y0)
            self.phi(:) = y0{1};
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Read and write to file
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function initializeStorage(self,group)
            group.addVariable(self.name,horzcat(self.model.wvt.spatialDimensionNames,'t'),isComplex=false,type="double",attributes=containers.Map({'isTracer'},{uint8(1)}));
        end

        function writeTimeStepToFile(self,group,outputIndex)
            group.variableWithName(self.name).setValueAlongDimensionAtIndex(self.phi,'t',outputIndex);
        end
    end

    methods (Static)
        function os = observingSystemFromGroup(group,model,outputGroup)
            %initialize a WVObservingSystem instance from NetCDF file
            %
            % Subclasses to should override this method to enable model
            % restarts. This method works in conjunction with -writeToFile
            % to provide restart capability.
            %
            % - Topic: Initialization
            % - Declaration: os = observingSystemFromGroup(group,wvt)
            % - Parameter model: the WVModel to be used
            % - Returns os: a new instance of WVObservingSystem
            arguments
                group NetCDFGroup {mustBeNonempty}
                model WVModel {mustBeNonempty}
                outputGroup WVModelOutputGroup
            end
            % most variables will be returned with this call, but we still
            % need to fetch (x,y,z), and the tracked variables
            vars = CAAnnotatedClass.requiredPropertiesFromGroup(group);

            parentGroup = outputGroup.group;
            nPoints = parentGroup.dimensionWithName("t").nPoints;
            vars.phi = parentGroup.readVariablesAtIndexAlongDimension('t',nPoints,vars.name);

            options = namedargs2cell(vars);
            os = WVTracer(model,options{:});
        end

        function vars = classRequiredPropertyNames()
            vars = {'name','isXYOnly','absTolerance','shouldAntialias'};
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            arguments (Output)
                propertyAnnotations CAPropertyAnnotation
            end
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
            propertyAnnotations(end+1) = CAPropertyAnnotation('name','name of Lagrangian particles');
            propertyAnnotations(end+1) = CANumericProperty('isXYOnly',{}, 'bool', 'whether the advection is only applied in x-y');
            propertyAnnotations(end+1) = CANumericProperty('absTolerance', {}, '','absolute tolerance of phi for the adaptive integrator');
            propertyAnnotations(end+1) = CANumericProperty('shouldAntialias',{}, 'bool', 'whether to antialias');
        end
    end
end