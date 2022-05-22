classdef ParticleFluxOperation < ModelOperation
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        name % ?
        outputVariables % Make private?
        f = []
    end

    methods
        function self = ParticleFluxOperation(outputVariables,f)
            self.outputVariables = outputVariables;
            for iVar=1:length(self.outputVariables)
                self.outputVariables(iVar).modelOp = self;
            end
            if nargin == 2
                self.f = f;
            end
        end

        function varargout = Compute(self,wvt,x,y,z)
            varargout = cell(1,length(self.outputVariables));
            [varargout{:}] = self.f(wvt,x,y,z);
        end
    end
end