classdef WVMooring < WVObservingSystem
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties (GetAccess=public, SetAccess=protected)
        isXYOnly
        phi
        absTolerance
        shouldAntialias
    end

    methods
        function self = WVMooring(model,options)
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
                options.nMoorings = 1
            end
            self@WVObservingSystem(model,"mooring");
            self.nFluxComponents = 0;

            N = options.nMoorings;

            % 1) "Ideal" continuous grid count in x
            fx = sqrt(N * model.wvt.Lx / model.wvt.Ly);

            % 2) Round to integer grid dims, ensure nx*ny ≥ N
            nx = max(1, round(fx));
            ny = max(1, ceil(N / nx));
            if nx*ny < N
                nx = max(1, ceil(fx));
                ny = max(1, ceil(N / nx));
            end
            for i=1:(nx*ny)

            end
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
            os = WVMooring(wvtX2,self.name);
            error('this needs to be implemented');
        end
    end
end