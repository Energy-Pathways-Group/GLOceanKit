classdef WVTracer < WVObservingSystem
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties (GetAccess=public, SetAccess=protected)
        isXYOnly
        phi
        absTolerance
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
                options.absTolerance = 1e-6; 
            end
            self@WVObservingSystem(model,options.name);
            if ~isfield(options.phi)
                error('You must specify the initial tracer field phi');
            end

            self.isXYOnly = options.isXYOnly;
            self.phi = options.phi;
            self.absTolerance = options.absTolerance;

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
                F = {-wvt.u .* wvt.diffX(self.phi) - wvt.v .* wvt.diffY(self.phi)};
            else
                F = {-wvt.u .* wvt.diffX(self.phi) - wvt.v .* wvt.diffY(self.phi) - wvt.w .* wvt.diffZF(self.phi)};
            end
        end

        function updateIntegratorValues(self,t,y0)
            self.phi = y0{1};
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Read and write to file
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function initializeStorage(self,group)
            group.addVariable(self.netCDFOutputTracers{iTracer},horzcat(self.model.wvt.spatialDimensionNames,'t'),type="double",attributes=containers.Map({'isTracer'},{true}));
        end

        function writeTimeStepToFile(self,group,outputIndex)
            group.variableWithName(self.name).setValueAlongDimensionAtIndex(self.phi,'t',outputIndex);
        end

        function os = observingSystemWithResolutionOfTransform(self,wvtX2)
            %create a new WVObservingSystem with a new resolution
            %
            % Subclasses to should override this method an implement the
            % correct logic.
            %
            % - Topic: Initialization
            % - Declaration: os = observingSystemWithResolutionOfTransform(self,wvtX2)
            % - Parameter wvtX2: the WVTransform with increased resolution
            % - Returns force: a new instance of WVObservingSystem
            os = WVTracer(wvtX2,self.name);
            error('this needs to be implemented');
        end
    end
end